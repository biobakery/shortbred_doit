#!/usr/bin/env python

import sys
import csv
import argparse
import re

import Bio
from Bio.Seq import Seq
from Bio import SeqIO


for seq in SeqIO.parse(sys.stdin, "fasta"):
	print seq.id, len(seq.seq)
    
    
