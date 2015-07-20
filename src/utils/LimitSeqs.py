#!/usr/bin/env python

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

import re
import sys
import datetime

import argparse


parser = argparse.ArgumentParser(description='Prints out first N sequences.')
parser.add_argument('-N', type=int, dest='iSeqs', help='Enter the number of sequences to print.', default=2000)
args = parser.parse_args()

def getSeqs(fileFasta):

    aSeqs = []

    for seq in SeqIO.parse(fileFasta, "fasta"):
        aSeqs.append(seq)

    return aSeqs


aSeqs = getSeqs(sys.stdin)

aSeqs = aSeqs[:args.iSeqs]

#print new list to fasta file
for seq in aSeqs:
    SeqIO.write( seq, sys.stdout, "fasta")
