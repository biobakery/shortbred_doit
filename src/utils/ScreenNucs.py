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

parser = argparse.ArgumentParser(description='Cut down nucleotide fasta file to sequences in corresponding protein file.')

parser.add_argument('--nucs', type=str, dest='sNucs', help='Enter the path and name of the nucleotide file.')
parser.add_argument('--prots', type=str, dest='sProts', help='Enter the path and name of the protein file.')
parser.add_argument('--log', type=str, dest='sLog', help='Enter the path of the log file.')
parser.add_argument('--out', type=str, dest='sOut', help='Enter the path of the output file.')

args = parser.parse_args()

def GetSeqNames(fileFasta):
	astrGeneNames = []
	for seq in SeqIO.parse(fileFasta, "fasta"):
		astrGeneNames.append(seq.id)
	return astrGeneNames

astrProts = GetSeqNames(args.sProts)
astrNucs = GetSeqNames(args.sNucs)

setProts = set(astrProts)
setNucs = set(astrNucs)

setIntersection = setProts.intersection(setNucs)

#load in nucs
dictNucs = {}
for seq in SeqIO.parse(args.sNucs, "fasta"):
	dictNucs[seq.id] = seq


fileLog = open(args.sLog,'w')
fileLog.write("Input Prots:")
fileLog.write(str(len(setProts)))
fileLog.write("Input Nucs:")
fileLog.write(str(len(setNucs)))
fileLog.write("Intersection:")
fileLog.write(str(len(setIntersection)))
fileLog.close()

ageneOut = []
for strName in setIntersection:
	ageneOut.append(dictNucs[strName])

SeqIO.write(ageneOut, args.sOut, "fasta")



