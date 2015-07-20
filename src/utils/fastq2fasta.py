#!/usr/bin/python
import sys

import Bio
from Bio.Seq import Seq
from Bio import SeqIO


for seq in SeqIO.parse(sys.stdin, "fastq"):
    SeqIO.write(seq, sys.stdout, "fasta")