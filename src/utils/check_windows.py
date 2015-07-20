#!/usr/bin/python

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

from collections import Counter

import re
import sys
import argparse


parser = argparse.ArgumentParser(description='Compares file of windows against orginial, clustered file to count windows, check for missing windows.')
parser.add_argument('--cf', type=file, dest='fClusteredFile', help='Enter the clustered file from the project.')
parser.add_argument('--list', type=bool,default=False, dest='bList', help='Set to \'True\' if you want the gene names.')
args = parser.parse_args()


def getSeqs(fileFasta):

    agSeqs = []

    for seq in SeqIO.parse(fileFasta, "fasta"):
        agSeqs.append(seq)

    return agSeqs

agClustSeqs = getSeqs(args.fClusteredFile)
agWindowSeqs = getSeqs(sys.stdin)

##############################################################################
#Prep Data

iTM = 0
iQM = 0
setHasTM = set()
setHasQM = set()

astrWindowNames = []
for i in agWindowSeqs:
    mtchName = re.search(r'(.*)(\_[TCQ]M[0-9]*)(\_\#[0-9]*)',i.id)
    strName = mtchName.group(1)
    astrWindowNames.append(strName)
    if mtchName.group(2)[:3] == "_TM":
        iTM +=1
        setHasTM.add(strName)
    elif mtchName.group(2)[:3] == "_QM":
        iQM +=1
        setHasQM.add(strName)


astrClustNames = []
for i in agClustSeqs:
    astrClustNames.append(i.id)

dictCounts = Counter(astrWindowNames)
atupCounts = dictCounts.items()
aCounts = []

for tup in atupCounts:
    aCounts.append(tup[1])

sClustNames = set(astrClustNames)
sWindowNames = set(astrWindowNames)

print "Number of True Markers:"
print iTM
print "Number of Quasi-Markers:"
print iQM
print len(agWindowSeqs)
print "Families with TM:"
print len(setHasTM)
print "Families with QM:"
print len(setHasQM)

print "Families missing markers (may have been dropped in quasi-clustering):"
print len(sClustNames.difference(sWindowNames))

dictCounts = Counter(aCounts)
print dictCounts.items()

print dictCounts[1]
print sum(dictCounts.values())-dictCounts[1]

if(args.bList):
    print "List of genes that have 0 windows:"
    for x in sClustNames.difference(sWindowNames):
        print x


"""
#######################################################
#Tabs on Windows file

print "Total number of markers:"
print len(agWindowSeqs)
print "True markers:"
print iTM
print "Quasi markers:"
print iQM

print "Breakdown of markers:"
print "(N, Number of genes with N markers)"
for x in Counter(aCounts).items():
    print x

#######################################################
#Check for missing genes



print "Unique genes in clustered file:"
print len(sClustNames)

print "Unique genes in markers file:"
print len(sWindowNames)


print "Number of genes in original clustered file that have 0 markers:"
print len(sClustNames.difference(sWindowNames))

if(args.bList):
    print "List of genes that have 0 windows:"
    for x in sClustNames.difference(sWindowNames):
        print x
"""

