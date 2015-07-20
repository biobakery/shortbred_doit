#! /usr/bin/env python

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

import sys, re, argparse
import math, random
from string import uppercase

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Data import CodonTable

from zopy.stats import weighted_chooser

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_strGoodAAs = re.sub( "[BOJUXZ]", "", uppercase )

c_dictSelfScore = {}
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

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------

def funcGetArgs ():
    """ """
    parser = argparse.ArgumentParser( formatter_class=argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument( "faa", help="fasta of protein sequences" )
    parser.add_argument( "--reads", default=1000, type=int, help="how many reads to make" )
    parser.add_argument( "--percentid", default=100, type=float, help="percent identity" )
    parser.add_argument( "--readlen", default=100, type=int, help="read length" )
    return parser.parse_args()

def funcLoadFAA( strPath, iFragLen ):
    """ load faa to dict; filter bad proteins """
    dictFAA = {}
    iBadCount = 0
    with open( strPath ) as fh:
        for record in SeqIO.parse( fh, "fasta" ):
            strSeq = str( record.seq )
            # exclude seqs matching a NOT good AA
            if not re.search( "[^%s]" % c_strGoodAAs, strSeq ) and len( strSeq ) >= iFragLen:
                dictFAA[record.name] = strSeq
            else:
                iBadCount += 1
    astrNames = dictFAA.keys()
    print >>sys.stderr, "skipped %d bad sequences (bad aa OR too short); %d remain" \
        % ( iBadCount, len( astrNames ) )
    return dictFAA

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

def funcReverseComplement ( strDNA ):
    """ reverse complement DNA sequence """
    return str( Seq( strDNA ).reverse_complement() )

def funcRandRead ( strProtSeq, strProtName, iRead, iReadLen, iFragLen, fMutRate ):
    """ create a single random read from a given protein sequence """
    iLen = len( strProtSeq )
    iStart = random.choice( range( iLen - iFragLen + 1 ) )
    strProtFrag = strProtSeq[iStart:iStart+iFragLen]
    # mutate?
    if fMutRate > 0:
        strProtFrag = funcMutateProtein( strProtFrag, fMutRate )
    strDNARead = funcProtToDNA( strProtFrag )
    iFrame = random.choice( range( 3 ) )
    strDNARead = strDNARead[iFrame:iFrame+iReadLen]
    # reverse complement?
    strStrand = "+"
    if random.random() < 0.5:
        strDNARead = funcReverseComplement( strDNARead )
        strStrand = "-"
    # format
    print ">protein_read|%08d|%04d|%s%d|%s" % ( iRead, iStart+1, strStrand, iFrame, strProtName )
    print strDNARead

# ---------------------------------------------------------------
# make reads
# ---------------------------------------------------------------

def funcMain ( ):
    args = funcGetArgs()
    # derived args
    strPathFAA = args.faa
    iTotalReads = args.reads
    fPercID = args.percentid
    if fPercID > 1:
        print >>sys.stderr, "interpret", fPercID, "as percent;",
        fPercID /= 100.0
        print >>sys.stderr, "fractional value is", fPercID
    fMutRate = 1 - fPercID
    iReadLen = args.readlen
    iFragLen = int( math.ceil( iReadLen / 3.0 ) )
    # load proteins
    dictFAA = funcLoadFAA( strPathFAA, iFragLen )
    dictProtLens = {k:len( v ) for k, v in dictFAA.items()}
    wc = weighted_chooser( dictProtLens )
    # make the reads (weighting by protein length)
    iRead = 0
    for strProtName in wc.iter_choice( iTotalReads ):
        iRead += 1
        strProtSeq = dictFAA[strProtName]
        funcRandRead( strProtSeq, strProtName, iRead, iReadLen, iFragLen, fMutRate )

if __name__ == "__main__":
    funcMain()
