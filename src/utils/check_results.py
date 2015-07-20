#!/usr/bin/env python

import sys
import csv
import argparse
import re

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

from collections import Counter

parser = argparse.ArgumentParser(description='Count matches in blast results and compare with wgs data to produce TP,FP,TN, and FN counts')
parser.add_argument('-fasta', default = "", type=str, dest='strFasta', help='Enter the path and name of your fasta file.')
parser.add_argument('-map', default = "", type=str, dest='strMap', help='Enter the path and name of your gene-centroid map.')



parser.add_argument('-endnum', default = "Y", type=str, dest='strEndNum', help='')
parser.add_argument('-filter', default = "N", type=str, dest='strFilter', help='Only count hits if they are above 90%')
args = parser.parse_args()

#This only works if each read a unique suffix!

#Get correct count from the fasta file

"""
Sample BLAST row:
VFG0676_238_408_0/1	VFG0676_TM_#01	93.94	33	2	0	1	99	15	47	1e-23	98.2	100
"""
    


dictWGS = {}
dictGeneMap = {}


iWGSReadCount = 0

#Make the gene-centroid map from uc file
for strLine in csv.reader(open(args.strMap),delimiter='\t'):
    if (strLine[0] == "H"):
        dictGeneMap[strLine[-2]] = strLine[-1]
    elif (strLine[0] == "C"):
        dictGeneMap[strLine[-2]] = strLine[-2]



#Load the WGS data. For each gene stub, copy in all the matching genes
for seq in SeqIO.parse(args.strFasta, "fasta"):
    mtchStub = re.search(r'(.*)(\_[0-9]*\_[0-9]*\_.*\/.*)',seq.id)
    dictWGS.setdefault(mtchStub.group(1),set()).add((seq.id))
    iWGSReadCount +=1
    


for key in dictWGS.keys():
    #Remove seqs that are not in the original prot file
    if key not in dictGeneMap.keys():
        #print key
        del dictWGS[key]
    #If gene!= cluster, add the entries for "gene" to the entries for "cluster"
    elif (key != dictGeneMap[key]):
        #print key, dictGeneMap[key]
        dictWGS[dictGeneMap[key]] =  dictWGS.setdefault(dictGeneMap[key],set()).union(dictWGS[key])
        del dictWGS[key]

#Read in the blast output: ("WGS","Sig/Seq")
dictBLAST = {}    
for aLine in csv.reader( sys.stdin, csv.excel_tab ):
    mtchStubWGS = re.search(r'(.*)(\_[0-9]*\_[0-9]*\_.*\/.*)',aLine[0])
    mtchStubSig = re.search(r'(.*)_.M[0-9]*_\#([0-9]*)',aLine[1])    
    if (args.strEndNum=="Y"):
        strStubSig = mtchStubSig.group(1)        
    else:
        strStubSig = aLine[1]
    
    if (args.strFilter=="Y"):
        fPct = float(aLine[2])
        iLength = int(aLine[3])
        if (fPct > 90.0):
            dictBLAST.setdefault(strStubSig,set()).add((aLine[0]))
        
            
    else:
        dictBLAST.setdefault(strStubSig,set()).add((aLine[0]))

#For the WGS data, group all the genes by their cluster

#dictGeneMap[key] = key's assigned cluster

"""
#print dictGeneMap["VFG1083"]
#print dictGeneMap["VFG0938"]
#print dictBLAST["VFG1083"]
#print dictWGS["VFG1083"]

#print len(dictWGS)
"""

"""
print "Length of dictWGS"
print len(dictWGS)
    


print dictBLAST["VFG1083"]
print dictWGS["VFG1083"]
"""   
     
astrResults = []
for strGene in dictWGS.keys():
    iTP = len(dictWGS[strGene].intersection(dictBLAST.get(strGene,set())))
    iFP = len(dictBLAST.get(strGene,set()).difference(dictWGS[strGene]))
    iFN = len(dictWGS[strGene].difference(dictBLAST.get(strGene,set())))
    iTN = iWGSReadCount - (iTP + iFP + iFN)
    iReal = iTP + iFN
    dAcc = float((iTP + iTN)) / iWGSReadCount
    
    
    dX = iFP/float(iFP + iTN)    
    dY = iTP/float(iTP + iFN)
    
    dAUC = dX*dY/2 + (1-dX)*dY + (1-dX)*(1-dY)/2
    
    astrResults.append([ strGene, str(iTP), str(iFP),  str(iTN),str(iFN),str(dAcc),str(dX),str(dY),str(dAUC)])
    
    #AUC = x*y/2 + (1-x)*y + (1-x)*(1-y)/2
    #= (lower left triangle) + (lower right rectangle) + (upper right triangle)

    #Here, x = FPR = FP/N, and y = TPR = TP/P

#print dictBLAST.items()
#print dictWGS["CAE55180"]
#print dictBLAST["CAE55180"]
astrHeader = ["Gene","TP","FP","TN","FN","Acc","FPR","TPR","AUC"] 
print "\t".join(astrHeader)
       
for strLine in astrResults:
    print "\t".join(strLine)
  

zipped = zip(*astrResults)
TPsum = sum(map(int,zipped[1]))
FPsum = sum(map(int,zipped[2]))
TNsum = sum(map(int,zipped[3]))
FNsum = sum(map(int,zipped[4]))


print "Total", TPsum, FPsum, TNsum, FNsum


