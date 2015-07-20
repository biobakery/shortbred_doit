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


parser.add_argument('--nucs', type=str, dest='sNucs', help='Enter the path and name of the nucleotide file.')
parser.add_argument('-N', type=int, dest='iN', help='Enter the number of genes to include in the file.',default =10)
parser.add_argument('--pct', type=float, dest='dPct', help='Enter the starting percentage, <= .30.', default=.20)
parser.add_argument('--log', type=str, dest='sLog',help='Enter the name of the log file')
parser.add_argument('--map', type=str, dest='sMap',help='Enter the name of the map file from usearch')


parser.add_argument('--props', type=str, dest='sProps', help='Enter the path and name of the file with gene names and proportions.')

args = parser.parse_args()


dictMap = {}
for aLine in csv.reader( open(args.sMap), csv.excel_tab ):
    dictMap[aLine[1]] = aLine[0]



atupPropData = []
dictNucs = {}
atupPrints = []


#Load in data from nuc file
for seq in SeqIO.parse(args.sNucs, "fasta"):
    dictNucs[seq.id] = seq

#Create an array of the percentages for the simulated genes

dPct = args.dPct
iN = args.iN
adProps = []

#First three percentages are handled here, rest are determined by loop.
adProps.append(dPct)
adProps.append(dPct/4)
adProps.append(dPct/4)

dRmndr = 1 - float(sum(adProps))

for i in range (3,iN-1):
    dRmndr = dRmndr/2
    adProps.append(dRmndr)
    
adProps.append(dRmndr)
adProps = sorted(adProps)


#Randomly select N genes from the nucs file

astrNames = random.sample(dictNucs.keys(), iN)


#Take the smallest percentage, that gets at least 10 prints of the gene


aaGenes = []

for i in range(0,iN):
    aGene = []
    aGene.append(astrNames[i])
    aGene.append(adProps[i])
    aaGenes.append(aGene)

iBase = 10    
dminPct = aaGenes[0][1]
iTotalGenes = 0
iSetLength = 0

for row in aaGenes:
    row.append(iBase * int((row[1]/dminPct)))
    iTotLength = row[2] * len(dictNucs[row[0]])    
    row.append(iTotLength)
    
    iTotalGenes += row[2]
    iSetLength += row[3]
    



for row in aaGenes:
    row.append(row[2]/float(iTotalGenes))
    row.append(row[3]/float(iSetLength))
    


fLog = open(args.sLog, 'w')

fLog.write( ("\t").join(["Clust Name","Gene Name", "User Pct", "Count", "LenXCount", "Rel Abundance (by Gene Count)", "Rel Abundance (by Gene Length)", "\n"]))




#print aaGenes[0]
#print iTotalGenes
#print iSetLength
   

for gene in aaGenes:
    iCount = gene[2]
    seq = dictNucs[gene[0]]
    for i in range(0, iCount):
        SeqIO.write(seq, sys.stdout, "fasta")
        
        
#Convert to string for Printing
for gene in aaGenes:
    gene[1] = str(round(gene[1],8))
    gene[2] = str(int(gene[2]))
    gene[3] = str(int(gene[3]))
    gene[4] = str(round(gene[4],8))
    gene[5] = str(round(gene[5],8))
    fLog.write(dictMap[gene[0]] + "\t")
    fLog.write( ("\t").join(gene))
    fLog.write("\n")

fLog.write("Total seqs: " + "\t"+ str(  iTotalGenes))
fLog.write("Total nucs: " + "\t"+ str(  iSetLength))


#want (genename, count)

"""
for tup in atupPrints:
    iCount = tup[1]
    while (iCount>0):
        SeqIO.write(dictNucs[tup[0]], sys.stdout, "fasta")
        iCount = iCount-1
"""
"""    
#Old code to load in data from props file
for astrLine in csv.reader( open(args.sProps), csv.excel_tab ):
    #aPropData.append(astrLine[0],astrLine[1])
    tupData = astrLine[0], astrLine[1]
    atupPropData.append(tupData)


atupPropData = sorted(atupPropData, key=lambda prop: prop[1])

strStart = atupPropData[0][0]
#print strStart, len(dictNucs[strStart])

iBase = 50* len(dictNucs[strStart])
dPctBase = float(atupPropData[0][1])

for row in atupPropData:
    strName = row[0]
    dPct =  float(row[1])
    iPrints = int(round(((dPct/dPctBase) * iBase)/ len(dictNucs[strName])))
    tupPrints = strName, iPrints
    atupPrints.append(tupPrints)

    
for tup in atupPrints:
    iCount = tup[1]
    while (iCount>0):
        SeqIO.write(dictNucs[tup[0]], sys.stdout, "fasta")
        iCount = iCount-1
"""


# Select 100 seqs randomly

#05 from Gene1
#10 from Gene2
#20 from Gene3
#40 from Gene4



#Look at the smallesr defined percentage 
    #Print that one 10 times




#Remaining 25 random from remaining 96 seqs


#Print gene1 50 times,
#print gene2 50*len(Gene1)/len(Gene1) times




