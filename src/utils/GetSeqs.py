import sys
import csv
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import re
import os


import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO


parser = argparse.ArgumentParser(description='This program pulls seqs with a given name out of a big fasta file .')




parser.add_argument('--name', type=str, dest='strName', help='Enter the name of the seq you want.')
parser.add_argument('--map', type=str, dest='strMap', help='Enter the name of the map you want to use.')
args = parser.parse_args()

dictFams = {}
for astrLine in csv.reader(open(args.strMap),delimiter='\t'):
	dictFams[astrLine[1]]=astrLine[0]

for seq in SeqIO.parse(sys.stdin, "fasta"):
	mtchGene = re.search(r'USR_(.*)_END',seq.id)
	if(mtchGene):
		strGeneName = str(mtchGene.group(1)).strip()
		#print strGeneName
	else:
		strGeneName =""

	if dictFams.get(strGeneName,"")==args.strName:
		SeqIO.write(seq, sys.stdout, "fasta")



