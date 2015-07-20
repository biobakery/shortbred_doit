#!/usr/bin/env python

import sys
import csv
import argparse
import re

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Drop sequences that are shorter or equal to N')
parser.add_argument('-N', default = 170, type=int, dest='iN', help='Cutoff length.')
args = parser.parse_args()


for seq in SeqIO.parse(sys.stdin, "fasta"):
    if (len(seq.seq) >= args.iN):
        SeqIO.write(seq, sys.stdout, "fasta")
    
    