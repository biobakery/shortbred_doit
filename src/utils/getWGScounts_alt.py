#!/usr/bin/env python

import re
import sys
import csv
import argparse
import os
import subprocess
import glob

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Count up user-added reads in wgs, group by prot family')
parser.add_argument('--map_in', type=str, dest='sMapIn',default="", help='Enter the path and name of the two column file connecting proteins to families.')
parser.add_argument('--wgs', type=str, dest='sWGS',default="", help='Enter the path and name of a file containing wgs.')

args = parser.parse_args()

dictFams={}
for astrLine in csv.reader( open(args.sMapIn,'r'), csv.excel_tab ):
    #Input file: family, protein
    #dictFams[protein]=[family]
	dictFams[str(astrLine[1]).strip()]=str(astrLine[0]).strip()



iReads = 0

dictGeneCounts = {}
dictFamCounts ={}

with open(args.sWGS,'r') as fCounts:
	fCounts.next()
	for astrLine in csv.reader( fCounts, csv.excel_tab ):
		mtch = re.search(r'USR_(.*)_END',astrLine[0])
		if(mtch):
			dictGeneCounts[mtch.group(1)]=float(astrLine[1])


#for key in dictFamCounts.keys():
#       print key + "\t" + str(dictFamCounts[key])


dictFinalCount = {}

for key in dictGeneCounts.keys():
	strFinalFam = dictFams.get(key,"NoCorrespondingProt")
	dictFinalCount[strFinalFam] = dictFinalCount.get(strFinalFam,0) + (dictGeneCounts[key])

print "Family" + "\t" + "Count"
for key in dictFinalCount.keys():
	print key + "\t" +  str(dictFinalCount[key])

#print "\n" + "There were " + str(iReads) + " user-supplied reads in the file."

#Debugging:
#print "Fams in dictionary: " , len(dictFams)
#print "Reads processed: ", iReads
