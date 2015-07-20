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

for strLine in (open(args.sWGS,'r')):
	mtch = re.search(r'USR_(.*)_END',strLine)
	#print mtch
	#print strLine
	#print mtch
    	if (mtch!=None):
		strID = mtch.group(1)
		#print strID

        	dictGeneCounts[strID] = dictGeneCounts.get(strID,0) + 1
		dictFamCounts[strID] = dictFamCounts.get(strID,0) + 1
        	iReads +=1

#for key in dictFamCounts.keys():
#       print key + "\t" + str(dictFamCounts[key])


for key in dictGeneCounts.keys():
	print key + "\t" + str(dictGeneCounts[key]) + "\t" + dictFams.get(key,"NoCorrespondingProt")

#print "\n" + "There were " + str(iReads) + " user-supplied reads in the file."

#Debugging:
#print "Fams in dictionary: " , len(dictFams)
#print "Reads processed: ", iReads
