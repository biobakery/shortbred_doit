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
parser.add_argument('--map', default = "", type=str, dest='sMap', help='Enter the path and name of your gene-centroid map.')
parser.add_argument('--wgs', default = "", type=str, dest='sWGS', help='Enter the path and name of your wgs file.')
parser.add_argument('--full', default = "", type=str, dest='sFull', help='Enter Y if using full sequences')


args = parser.parse_args()

#This only works if each read a unique suffix!

#Get correct count from the fasta file

"""
Sample BLAST row:
VFG0676_238_408_0/1	VFG0676_TM_#01	93.94	33	2	0	1	99	15	47	1e-23	98.2	100
"""
    



dictGeneMap = {}
for aLine in csv.reader( open(args.sMap), csv.excel_tab ):
    dictGeneMap[aLine[1]] = aLine[0]



dictPositives = {} 
dictWGSCounts = {}

iTP = 0
iFP = 0

for key in dictGeneMap.values():
    dictPositives[key] = [iTP,iFP]    

#Make a dict for wgs genes in our set of interest
#Count EVERY read for the grand total

iWGSReadCount = 0

for seq in SeqIO.parse(open(args.sWGS), "fasta"):
    mtchStubWGS = re.search(r'(.*)(\_[0-9]*\_[0-9]*\_.*\/.*)',seq.id)
    strGeneID = mtchStubWGS.group(1)
    if strGeneID in dictGeneMap.keys():
        strFamID = dictGeneMap[strGeneID]
        dictWGSCounts[strFamID] = dictWGSCounts.get(strFamID,0) + 1
    
    iWGSReadCount += 1

for aLine in csv.reader( sys.stdin, csv.excel_tab ):
    
    mtchStubWGS = re.search(r'(.*)(\_[0-9]*\_[0-9]*\_.*\/.*)',aLine[0])
    strStubWGS = mtchStubWGS.group(1)        
    
    if (args.sFull == "Y"):
        strStubSig = aLine[1]
    else:
        mtchStubSig = re.search(r'(.*)_.M[0-9]*_\#([0-9]*)',aLine[1])            
        strStubSig = mtchStubSig.group(1)        
    
    
    strSigFam = dictGeneMap.get(strStubSig,"NoSigFam")    
    strWGSFam = dictGeneMap.get(strStubWGS,"NoWGSFam")
    
    if (strSigFam==strWGSFam):
        dictPositives[strSigFam] = [dictPositives[strSigFam][0]+1,dictPositives[strSigFam][1] ]
    else:
        dictPositives[strSigFam] = [dictPositives[strSigFam][0],dictPositives[strSigFam][1]+1 ]


print "Family" + "\t" + "TP" +"\t" + "FP" + "\t" + "TN" + "\t" + "FN"

for key in dictPositives.keys():
    iTotal = dictWGSCounts.get(key,0)    
    iTP = dictPositives[key][0]
    iFP = dictPositives[key][1]
    
    iFN = iTotal - iTP
    iTN = iWGSReadCount - iFN - dictPositives[key][0] - dictPositives[key][1]
    
    print key + "\t" + str(iTP) +"\t" + str(iFP) + "\t" + str(iTN) + "\t" + str(iFN)        
        

#print dictPositives["VFG1368"] 


#Checking
   
    #If both stubs have the sam fam, add one to the TP

    #Else, add one to the FP


