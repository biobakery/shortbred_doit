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


parser = argparse.ArgumentParser(description='This program compares two fasta files and reports how they differ .')


parser.add_argument('--f1', type=str, dest='strF1', help='Fasta1.')
parser.add_argument('--f2', type=str, dest='strF2', help='Fasta2.')
args = parser.parse_args()

dictF1 = {}
dictF2 = {}

for seq in SeqIO.parse(args.strF1, "fasta"):
    dictF1[seq.id] = seq
for seq in SeqIO.parse(args.strF2, "fasta"):
    dictF2[seq.id] = seq

setF1 = set(dictF1.keys())
setF2 = set(dictF2.keys())
print "\n\nSequences in " + args.strF1 + " not in " + args.strF2 + "\n"
print(setF1.difference(setF2))
print "\n\nSequences in " + args.strF2 + " not in " + args.strF1 + "\n"
print(setF2.difference(setF1))

print "\n" + args.strF1 + "\n" + "*****************"
for strID in setF1.difference(setF2):
    SeqIO.write(dictF1[strID], sys.stdout, "fasta")

print args.strF2 + "\n" + "*****************"
for strID in setF2.difference(setF1):
    SeqIO.write(dictF2[strID], sys.stdout, "fasta")

print "Check sequences in both files. Print if their material differs."
for strID in setF1.intersection(setF2):
    if dictF1[strID].seq!=dictF2[strID].seq:
        print "File1:"
        SeqIO.write(dictF1[strID], sys.stdout, "fasta")
        print "File2:"
        SeqIO.write(dictF2[strID], sys.stdout, "fasta")

