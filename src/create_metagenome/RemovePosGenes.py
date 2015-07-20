#!/usr/bin/env python

#Jim Kaminski
#Huttenhower Lab
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

import sys
import csv
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import re
import os
import datetime

import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

c_strUSEARCH	= "/n/home11/jkaminski/usearch/usearch"
c_strBLASTX	= "/n/home11/jkaminski/blast/ncbi-blast-2.2.28+/bin/blastx"


##############################################################################
def check_create_dir( strDir ):

	if not os.path.exists( strDir ):
		os.makedirs( strDir )
	return strDir

def RunUSEARCH (strWGS,strBlastOut, strDB):

	subprocess.check_call([c_strUSEARCH, "--usearch_local", strWGS, "--db", strDB,
		"--id", str(args.dID),"--userout", strBlastOut,
		"--threads", str(args.iThreads), "--userfields", "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl+qs+ts"])
	return

def getOverlapCountsMOD (fileBlast, dIDcutoff, iRegionLength, dictOverlap,iOffset=0):
#Makes a dictionary of form (GeneID, [0,0,0,0,1,1...]), where the number indicates
#the number of times an amino acid overlaps with a region in the blast output.

#Read in the blast output line by line
#If line meets criteria, update overlap region

    iGeneCount = 0
    iLine =0

    for aLine in csv.reader( open(fileBlast, 'rb'), delimiter='\t' ):

        iLine+=1

        strQueryID = aLine[0]
        strSubId = aLine[1]
        dIdentity =float(aLine[2])/100
        iAln= int(aLine[3] )
        iMismatch = int(aLine[4])
        iGap = int(aLine[5] )
        iQStart = int(aLine[6] )
        iQEnd = int(aLine[7])
        iSubStart = int(aLine[8] )
        iSubEnd= int(aLine[9] )
        deVal = float(aLine[10] )
        dBit = float(aLine[11])
        iQLength = int(aLine[12])


        if (dIdentity >= dIDcutoff) and (iAln >= iRegionLength):
            for i in range(iQStart-1+iOffset, iQEnd-iOffset):
                dictOverlap[strQueryID][i]=dictOverlap[strQueryID][i]+1

	"""
    #Once the loop is done, remember to add the last window.
    dictAAOverlapCounts.setdefault(strCurQuery.strip(), aiCounts)
    iGeneCount = iGeneCount+1
	"""
    return dictOverlap

#################################################################################
parser = argparse.ArgumentParser(description='Remove Positive Genes  \n \
This program was built to examine the bacterial genomes used in ShortBRED\'s evaluation \
and remove ARDB and VF protein sequences before putting them into synthetic genomes.',formatter_class=RawTextHelpFormatter)

parser.add_argument('--prots', type=str, dest='strProts',default= "", help='Enter the path and name of the proteins you want to remove from the genome \
typically the ARDB or VF prot fasta file.')
parser.add_argument('--tmpgenome', type=str, dest='strTmpGenome',default= "", help='Make a copy of the genome with names BioPython can handle.')
parser.add_argument('--genome', type=str, dest='strGenome',default= "", help='Enter the path and name of the bacterial genome fasta file.')
parser.add_argument('--clean', type=str, dest='strCleanOut',default= "", help='Enter the path and name to output the clean fasta file.')

parser.add_argument('--id', type=float, dest='dID', help='Enter the percent identity for USEARCH', default = .95)
parser.add_argument('--threads', type=int, dest='iThreads', help='Enter the number of CPUs available for USEARCH.', default=3)
parser.add_argument('--len', default = .10, type=float, dest='dL', help='Enter the length maximum for IDing a region. l=(length hit region)/(length query gene) Examples: .30, .20, .10,... ')
parser.add_argument('--tmp', default = "tmp", type=str, dest='sTmp', help='Enter the name of the tmp directory')

args = parser.parse_args()


dirTmp = args.sTmp
dirTmp = check_create_dir( dirTmp)


# These next few lines make the KEGG genome readable by BioPython and USEARCH.
# USEARCH problem:
# USEARCH can't handle very long genes, so genes that are longer than iMaxGeneLength
# will get broken up into smaller genes.

# BioPython problem:
# BioPython was only getting the first set of letters from the gene header, which
# was causing it to think it was reading the same gene twice.

# Example of problem:
#	KEGG Name: ">gn:cdf [NC_009089] Clostridium difficile 630 chromosome, complete genome."
#	BioPython translation ">gn:cdf" --> causes repeats in the dictionary

streamGenome = open(args.strGenome,'r')
streamTmpGenome = open(args.strTmpGenome,'w')
iMaxGeneLineLength = 800
iLineCount = 0
iGeneCount = 0
strName = ""

for strLine in streamGenome:
	iLineCount+=1
	if strLine.find(">")==0:
		iGeneCount+=1
		iLineCount=0
		strLine = strLine[:7] + "_" + str(iGeneCount).zfill(3) + '\n'
		strName = strLine[1:7]

	if iLineCount > iMaxGeneLineLength:
		iGeneCount+=1
		iLineCount = 0
		strLine = ">" + strName + "_brk_" + str(iGeneCount).zfill(3) + '\n'
	streamTmpGenome.write(strLine)
streamTmpGenome.close()
streamGenome.close()

#Make overlap windows
dictOverlap = {}
dictGenes = {}
for gene in SeqIO.parse(args.strTmpGenome, "fasta"):
	dictOverlap[gene.id] = len(gene)*[0]
	dictGenes[gene.id] = gene



#Make usearch db
strDBName = str(dirTmp) + os.sep + os.path.basename(str(args.strProts)) + ".udb"
p = subprocess.check_call([c_strUSEARCH, "--makeudb_usearch", args.strProts,
	"--output", strDBName])

#Call USEARCH
strBlastOut = dirTmp + os.sep + "hits.out"
RunUSEARCH (args.strTmpGenome,strBlastOut, strDBName)

#Find overlap
dictOverlap = getOverlapCountsMOD(fileBlast = strBlastOut, dIDcutoff = args.dID, iRegionLength=0, dictOverlap=dictOverlap,iOffset=0)

#Replace bacterial nuc's that are really positives (ARDB or VF) with "N"'s
for key in dictGenes.keys():
	aiWindow = dictOverlap[key]
	astrGene = list(dictGenes[key])

	for i in range(len(aiWindow)):
		if aiWindow[i] >= 1:
			astrGene[i] = "N"
	dictGenes[key] = SeqRecord(Seq("".join(astrGene)),id = key, description = "")

streamOut = open(args.strCleanOut,'w')

for key in dictGenes.keys():
	SeqIO.write(dictGenes[key], streamOut, "fasta")

streamOut.close()
