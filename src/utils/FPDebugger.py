#!/usr/bin/env python
#####################################################################################
#Copyright (C) <2013> Jim Kaminski and the Huttenhower Lab
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# This file is a component of ShortBRED (Short, Better REad Database)
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Jim Kaminski, jjk451@mail.harvard.edu).
#####################################################################################

#Jim Kaminski
#Huttenhower Lab


import sys
import csv
import argparse
from argparse import RawTextHelpFormatter
import re
import os
import bisect

import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO



# Load in Markers, save in dict
# Load in Input Prots, save in dict
# Load in family dictionary

# Read in text file of Marker FP's
# For each pair of (Read, Marker):
#   Print what family the original input seq belongs to
#	Print marker seq
#   Print original input seq

#   Print *ALL* hits from appropriate BLAST txt (obtain from ShortBRED-Identify's tmp folder)

parser = argparse.ArgumentParser(description='FPDebugger - This program outputs the prot sequences, associated families, \
and BLAST hits from ShortBRED-Identify for FP WGS_Read-Marker matches',formatter_class=RawTextHelpFormatter)


parser.add_argument('--FP', type=str, dest='strFP', help='Enter the path and name of the FP file.')
parser.add_argument('--markers', type=str, dest='strMarkers', help='Enter the path and name of the marker file (fasta).')
parser.add_argument('--prots', type=str, dest='strProts', help='Enter the path and name of the prots file (fasta).')
parser.add_argument('--fams', type=str, dest='strFams', help='Enter the path and name of the family file.')
parser.add_argument('--blast', type=str, dest='strBlast', help='Enter the path and name of the ShortBRED-Identify BLAST data (selfblast or refblast).')
parser.add_argument('--table', type=str, dest='strTable', help='Enter the path and name of the Marker FP table.')

args = parser.parse_args()


##############################################################################
def PrintBlastHits(aastrBlastResults, strMarker, strRead,seqMarker):
	csvwResults = csv.writer(sys.stdout,delimiter="\t")
	#csvwResults.writerow(["MarkerFam","ProtFam","ID","AlnLength","Mismatches","Gap","MFam_Start","MFam_End","RFam_Start","MFam_End"])
	csvwResults.writerow(["Hits and Section of Original Query Prot"])
	csvwResults.writerow(["MarkerFam","ProtFam","ID","PctLen"])

	#Get the first column of the sorted results (the Query prot's name)
	astrIndex = zip(*aastrBlastResults)[0]

	if strMarker in astrIndex:
		iRow = index(astrIndex, strMarker)
		print iRow, len(aastrBlastResults), strMarker
		bSearching = True


		while(bSearching == True and iRow < len (aastrBlastResults)):
			aLine = aastrBlastResults[iRow]
			if aLine[0]==strMarker and aLine[1]==strRead:
		 		iQStart = int(aLine[6] )
				iQEnd   = int(aLine[7] )
				csvwResults.writerow( [aLine[0],aLine[1],aLine[2],"%0.2f" % (int(aLine[3])/float(len(seqMarker)))])
				print seqMarker.seq[iQStart:iQEnd]
				iRow +=1
			elif (aastrBlastResults[iRow][0]==strMarker and aLine[1]!=strRead):
				iRow +=1
			else:
				bSearching = False
	else:
		print "NO BLAST RESULTS FOR PROTEIN"

        #csvwResults.writerow

	"""
	for aLine in csv.reader( open(fileBlast, 'rb'), delimiter='\t' ):
		iQStart = int(aLine[6] )
		iQEnd   = int(aLine[7] )
		if aLine[0]==strMarker and aLine[1]==strRead:
			csvwResults.writerow( [aLine[0],aLine[1],aLine[2],"%0.2f" % (int(aLine[3])/float(len(seqMarker)))])
			print seqMarker.seq[iQStart:iQEnd]
            #csvwResults.writerow
	"""
###############################################################################
#COPIED and MODIFIED FROM http://docs.python.org/2/library/bisect.html
def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
	print "Problem with ", x
    raise ValueError

################################################################################
#Load Markers (fasta), InputProts (fasta), and families (two column txt file[fam,seq])

dictMarkers = {}
dictMarkerFPs = {}
for seq in SeqIO.parse(args.strMarkers, "fasta"):
		dictMarkers[seq.id]= seq
		dictMarkerFPs[seq.id] = 0

dictProts   = {}
for seq in SeqIO.parse(args.strProts, "fasta"):
        dictProts[seq.id]= seq

dictFams    = {}
for astrLine in csv.reader(open(args.strFams),delimiter='\t'):
	dictFams[astrLine[1]]=astrLine[0]

#Load Blast Results, sort by first column
aastrBlast = []
for aLine in csv.reader( open(args.strBlast, 'rb'), delimiter='\t' ):
	aastrBlast.append(aLine)

aastrBlast.sort(key=lambda col: col[0])




##################################################################################
# Read in text file of Marker FP's
# For each pair of (Read, Marker):
#   Print what family the original input seq belongs to
#	Print marker seq
#   Print original input seq
iHitCount = 0
for aLine in csv.reader( open(args.strFP, 'rb'), delimiter='\t' ):


	strRead = aLine[0]
	strMarker = aLine[1]
	dictMarkerFPs[strMarker] = dictMarkerFPs[strMarker]+1


	mtchProt = re.search(r'USR_(.*)_END',strRead)
	if(mtchProt):
		iHitCount +=1
		strReadProt = mtchProt.group(1).strip()

		mtchMarkerFam = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',strMarker)
		strMarkerFam = mtchMarkerFam.group(1).strip()

		print "**********************************************************"
		print "Bad match: " + strMarkerFam +" to " + strReadProt  + ", #" +str(iHitCount).zfill(4)

		print "Marker:", strMarker
		print "Marker Family:", strMarkerFam
		print "Read:", strReadProt
		print "Read Family:", dictFams[strReadProt], '\n','\n'

		SeqIO.write(dictMarkers[strMarker], sys.stdout,"fasta")
		print '\n'
		SeqIO.write(dictProts[strReadProt], sys.stdout,"fasta")
	 	print '\n'

		PrintBlastHits(aastrBlast, strMarkerFam, dictFams[strReadProt],dictProts[strMarkerFam])
	else:
		print "**********************************************************"
		print "Bad match: Bacterial FP"
		iHitCount +=1

# atupMarkerFPs = [ (AAA33_TM_#01, 23),...]
atupMarkerFPs = dictMarkerFPs.items()
atupMarkerFPs.sort(key=lambda x: x[1],reverse=True)

#Print out table of Marker FP's
with open(args.strTable, 'w') as csvfileMarker:
	csvwMarkerResults = csv.writer( csvfileMarker, csv.excel_tab )
	csvwMarkerResults.writerow(["Marker","FPs", "Marker Length"])
	for tup in atupMarkerFPs:
		if tup[1]>0:
			csvwMarkerResults.writerow([tup[0],tup[1],len(dictMarkers[tup[0]])*3])
