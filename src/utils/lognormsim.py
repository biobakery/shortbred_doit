#!/usr/bin/env python
"""
Created on Mon Oct  8 10:13:40 2012

@author: jim
"""

import sys
import argparse
import subprocess
import csv
import re
import random
import os

import Bio
from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers, and produces a relative abundance table.')

#INPUTS
parser.add_argument('--nucs', type=str, dest='sNucs', help='Enter the path and name of the input genes.')
parser.add_argument('--ctrlnucs', type=str, dest='sCtrlNucs', help='Enter the path and name of the control genes.')
parser.add_argument('--map', type=str, dest='sMap',help='Enter the name of the map file connecting genes to centroids')

#PARAMETERS
parser.add_argument('-N', type=int, dest='iN', help='Enter the number of input genes to include in the file.',default =10)
parser.add_argument('--CN', type=int, dest='iCN', help='Enter the number of control genes to include in the file.',default =10)
parser.add_argument('--min', type=int, dest='iMinLength', help='Enter the minimum length in bp for a gene to be used.',default =10)

#OUTPUT
#(fasta goes to stdout)
parser.add_argument('--log', type=str,dest='sLog',default="log.txt",help='Enter the name of the log file')
parser.add_argument('--byclust', type=str, dest='sFamLog', default="logclust.txt",help='Enter name for output file, props summed by cluster family')

args = parser.parse_args()

#######################################################################
#FUNCTIONS
########################################################################
def PullSeqs (strFileName, iN,setIn = set()):
    dictNucs = {}    
    #Load in data from nuc file
    for seq in SeqIO.parse(strFileName, "fasta"):
        dictNucs[seq.id] = seq
        
    
    #Strip out  seqs shorter than minLength
    for key in dictNucs.keys():
        if (len(dictNucs[key]) < args.iMinLength):
            del dictNucs[key]
    
    if (len(setIn)>0):
        for key in dictNucs.keys():
            if key not in setIn:
                del dictNucs[key]
    
    
    #Randomly select N genes from the nucs file
    astrNames = random.sample(dictNucs.keys(), iN)

    dictKeep = {}
    for strID in astrNames:
        dictKeep[strID] = dictNucs[strID]
        
    return dictKeep

################################################################
#PROGRAM
################################################################

#Get dict of (prot,cluster_name)
dictMap = {}
for aLine in csv.reader( open(args.sMap), csv.excel_tab ):
    dictMap[aLine[1]] = aLine[0]

setProts = set(dictMap.keys())

#Get dicts of (seqname, seq)
dictNucs = PullSeqs(args.sNucs,args.iN,setProts)
dictCtrls = PullSeqs(args.sCtrlNucs,args.iCN)



#Mean for lognomral distribution, sigma = 1/2 * Mu
iMu = 2
#iSigma = 1

dictClustCounts = {}
aaGeneCounts = []
iSum = 0
for key in dictNucs.keys():
    iCount = int(random.lognormvariate(iMu,iMu/2 ))
    aRow = [key, str(iCount)]
    iSum += iCount

    #Sum by protein family, add to dict Clust Counts    
    strFam = dictMap[key]
    dictClustCounts[strFam] = dictClustCounts.get(strFam,0) + iCount    
    
    aaGeneCounts.append(aRow)
    for i in range(0,iCount):    
        SeqIO.write(dictNucs[key], sys.stdout, "fasta")        

    
        
for key in dictCtrls.keys():
    iCount = int(random.lognormvariate(iMu,iMu/2 ))
    aRow = [key, iCount]
    iSum += iCount
    
    aaGeneCounts.append(aRow)
    for i in range(0,iCount):    
        SeqIO.write( dictCtrls[key], sys.stdout, "fasta")  



f = open(args.sLog, 'w')

for aRow in aaGeneCounts:
    strName = str(aRow[0])
    dProp = float(float(aRow[1])/iSum)
    f.write(strName +'\t' + "{0:.6f}".format(dProp) + '\n')
    
f.close()

    
f = open(args.sFamLog, 'w')

for key in dictClustCounts:
    dProp = float(dictClustCounts[key])/iSum
    f.write(key +'\t' + "{0:.6f}".format(dProp) + '\n')

f.close()
    

