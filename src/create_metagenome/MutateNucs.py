#!/usr/bin/env python
"""
Jim Kaminski - 6/29/2015
This code is a modified version of  Eric Franzosa's scriptProteinReads.py,  
whic can be found on bitbucket at https://bitbucket.org/franzosa/zopy/. 
The original is located in scripts/scriptProteinReads.py .

Some of the ARDB proteins we use have an "X" in them, so I set their 
self score to be 0. Since they are ambiguous, I think it is fair to 
change them to anything. I also assigned it a random choice from the
existing codons.

"""
"""
Make DNA reads from proteins
===============================================
* Loads proteins from FASTA
* Samples peptides (~reads) proportional to length
* Selects sites to mutate inversely proportional to exp( blosum62[x][x] )
* Chooses new residue proportional to exp( blosum62[x][y] ); x != y
* Converts mutated peptide to compatible DNA sequence
* Picks random reading frame
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import sys, re, argparse, os
#sys.stderr.write("MutateNucs.py is using the following for its path:")
#sys.stderr.write(" ".join(os.environ['PYTHONPATH'].split(os.pathsep)))
import math, random
from string import uppercase

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Data import CodonTable
from Bio.Alphabet import generic_dna

from zopy.stats import weighted_chooser


parser = argparse.ArgumentParser(description='Mutate Protein Sequences \n This program takes amino acid sequences, mutates them, and generates corresponding nucleotide sequences.')


parser.add_argument('--mutrate', default = .03, type=float, dest='dMutRate', help='Enter the pct of AA\'s to mutate. Examples: .03, .10, ,... ')

args = parser.parse_args()


# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------
 
c_strGoodAAs = re.sub( "[BOJUXZ]", "", uppercase )

c_dictSelfScore = {}
c_dictSelfScore["X"]  = 0
c_dictdictExchangeScores = {}

for ( strAA1, strAA2 ), fValue in blosum62.items():
    if strAA1 in c_strGoodAAs and strAA2 in c_strGoodAAs:
        if strAA1 != strAA2:
            c_dictdictExchangeScores.setdefault( strAA1, {} )[strAA2] = math.exp( fValue )
            c_dictdictExchangeScores.setdefault( strAA2, {} )[strAA1] = math.exp( fValue )
        else:
            # higher self score should be mutated less, hence negative
            c_dictSelfScore[strAA1] = math.exp( -fValue )
# convert exchanges to weighted chooser objects
c_dictChoosers = {}
for strAA, dictWeights in c_dictdictExchangeScores.items():
    c_dictChoosers[strAA] = weighted_chooser( dictWeights )

# process codon table
c_dictCodons = {}
for strCodon, strAA in \
        CodonTable.unambiguous_dna_by_name["Standard"].forward_table.items():
    if strAA in c_strGoodAAs:
        c_dictCodons.setdefault( strAA, [] ).append( strCodon )


print 

astrRandCodons = c_dictCodons[random.choice(c_dictCodons.keys())]
c_dictCodons["X"] = random.choice(astrRandCodons)

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------

def funcMutateProtein ( strSeq, fMutRate ):
    """ mutate the sequence string """
    iLen = len( strSeq )
    iMuts = int( fMutRate * iLen )
    astrSeq = [k for k in strSeq]
    # part 1: pick the sites to mutate based on their relative mutability
    dictSiteWeights = {k:c_dictSelfScore[astrSeq[k]] for k in range( iLen )}
    wc = weighted_chooser( dictSiteWeights )
    dictSitesToMutate = {}
    while len( dictSitesToMutate ) < iMuts:
        dictSitesToMutate[wc.choice()] = 1
    # part 2: make sensible mutations at those sites
    for iSite in dictSitesToMutate:
        astrSeq[iSite] = c_dictChoosers[astrSeq[iSite]].choice()
    strNewSeq = "".join( astrSeq )
    return strNewSeq



def funcProtToDNA ( strProtSeq ):
    """ create non-biological CDS by joining random compatible codons """
    return "".join( [random.choice( c_dictCodons[strAA] ) for strAA in strProtSeq] )
    
    
# ---------------------------------------------------------------
# Mutate input sequences and print results.
# ---------------------------------------------------------------
    
for seqProt in SeqIO.parse(sys.stdin, "fasta") :
        seqProt.seq = Seq(funcMutateProtein(seqProt.seq.upper(), args.dMutRate ),seqProt.seq.alphabet )
        seqNuc = seqProt
        seqNuc.seq = Seq(funcProtToDNA(seqProt.seq),generic_dna)
        SeqIO.write(seqNuc, sys.stdout, "fasta")
