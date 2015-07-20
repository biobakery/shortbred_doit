#!/usr/bin/env python

import sys
import argparse
import subprocess
import csv
import re
import os

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers, and produces a relative abundance table.')


parser.add_argument('--fasta', type=str, dest='sFasta', help='Enter the path and name of the genes of interest file (proteins).')
parser.add_argument('--map', type=str, dest='sMap', help='Enter the path and name of the map.')


args = parser.parse_args()




#Load gene-family map
dictGeneMap = {}
for aLine in csv.reader( open(args.sMap), csv.excel_tab ):
    dictGeneMap[aLine[1]] = aLine[0]

#Load the fasta file, get seq length (used for normalizing results)
dictMarkerLen = {}
for seq in SeqIO.parse(args.sFasta, "fasta"):
	strStub = seq.id
	dictMarkerLen[strStub] = len(seq) 
    
#Go through the blast hits, sum hits by family
dictBLAST = {}    
for aLine in csv.reader( sys.stdin, csv.excel_tab ):
    strProtFamily = aLine[1]
    dictBLAST[strProtFamily] = dictBLAST.setdefault(strProtFamily,0)+1
	

for strProt in dictBLAST.keys():
    print strProt + "\t" +  str(float(dictBLAST[strProt])/dictMarkerLen[strProt])
    
