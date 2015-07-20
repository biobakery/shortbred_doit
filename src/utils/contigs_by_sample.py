#!/usr/bin/env python

import sys
import csv
import argparse
import re

import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='This program takes the fna file of contigs and produces a set of \
organizes them into fasta files for each sample + resistance combination .')

parser.add_argument('--contigs', type=str, dest='strName', help='Big file of contigs.')
parser.add_argument('--dir', type=str, dest='strDir', help='Directory for output.')
args = parser.parse_args()

iCount = 0
iTotalBases = 0
dictCounts  = {}

for seq in SeqIO.parse(args.strName, "fasta"):
	iCount+= 1
	iTotalBases+= len(seq)
	astrName = re.split("_",seq.id)
	strName = astrName[0] +"_" + astrName[1]
	#print strName
	strFile = args.strDir + strName + ".fna"
	if strName in dictCounts:
		dictCounts[strName] = dictCounts[strName]+1
		with open(strFile,'a') as fileContigs:
			#seq = SeqRecord(seq.seq.translate(stop_symbol="X"),id=seq.id)
			SeqIO.write(seq, fileContigs,"fasta")
	else:
		dictCounts[strName]= 1
		with open(strFile,'w') as fileContigs:
					#seq = SeqRecord(seq.seq.translate(stop_symbol="X"),id=seq.id)
					SeqIO.write(seq, fileContigs,"fasta")


for strName in sorted(dictCounts.keys()):
	print strName, dictCounts[strName]

print "Average length: "
print str(iTotalBases/float(iCount))

