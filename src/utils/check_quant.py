#!/usr/bin/env python

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


parser = argparse.ArgumentParser(description='This program merges results from ShortBRED quantify and simbase.csv .')

parser.add_argument('--simbase', type=str, dest='sSim', help='Enter the path and name of the simulation proportions file.')
parser.add_argument('--results', type=str, dest='sResults', help='Enter the path and name of the results file.')
args = parser.parse_args()






aaSimBase = []
setSimNames = set()
iLine=0

#open simbase file, drop the header
for aLine in csv.reader( open(args.sSim), csv.excel_tab ):
    if (iLine!=0 and len(aLine)>5):
        aaSimBase.append(aLine)
        setSimNames.add(aLine[0])
    iLine+=1

#Load Results
aaResults = []




for aLine in csv.reader( open(args.sResults), csv.excel_tab ):
    if (iLine!=0 and len(aLine)==2 and aLine[0]!="" and aLine[1]!=""):
        aaResults.append(aLine)

print aaResults    
dResTotal = sum([float(row[1]) for row in aaResults])
    

print ("\t").join(["Clust Name", "Rel Abundance (by Gene Count)", "Rel Abundance (by Gene Length)", "ShortBRED"])

#combine data for clusters    
for name in sorted(setSimNames):
    
    #aData = [row for row in aaSimBase if row[0]==name]
    dRelCount = sum([float(row[5]) for row in aaSimBase if row[0]==name])
    dRelLength =sum([float(row[6]) for row in aaSimBase if row[0]==name])
    
    dSBValue = sum([float(row[1]) for row in aaResults if row[0]==name])
    
    strResults = [name, str(dRelCount), str(dRelLength), str(dSBValue/dResTotal)]
    print ("\t").join(strResults)

#print aaResults

aaNotInSim = []

for row in aaResults:
    if row[0] not in setSimNames:
        aaNotInSim.append(row)
        
dTotOther = sum([float(row[1]) for row in aaNotInSim])
strResults =["Other", "", "", str(dTotOther/dResTotal)]
print ("\t").join(strResults)

print "\n", "\n"

#Print results
for aLine in open(args.sResults):
    print aLine
print "\n", "\n"
#Print simbase.
for aLine in  open(args.sSim):
    print aLine