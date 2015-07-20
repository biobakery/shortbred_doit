#! /usr/bin/env python

import os, sys, re, glob, argparse
from numpy import median
from zopy.utils import reader

# argument parsing ( python argparse )
parser = argparse.ArgumentParser()
parser.add_argument( '-i', '--input', default=None, help='filename [stdin if empty]' )
parser.add_argument( '-k', '--ksize', default=5, type=int, help='kmer size [default 5]' )
parser.add_argument( '-w', '--words', action="store_true", help='use words instead of kmers' )
parser.add_argument( '-c', '--mincount', default=10, type=int, help='min kmer count [default 10]' )
parser.add_argument( '-e', '--extent', default=0.1, type=float, help='how deep to got into ranked list [default 10% from each end]' )
parser.add_argument( '-j', '--minmerge', default=0.9, help='local jaccard similarity for cluster merger' )
parser.add_argument( '-s', '--simplify', action="store_true", help='lowercase, remove all but [a-z0-9]' )
args = parser.parse_args()

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

# constants from above
mincount = args.mincount
k        = args.ksize
extent   = args.extent
minmerge = args.minmerge
simplify = args.simplify

# constants
fh = open( args.input ) if args.input is not None else sys.stdin

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def overlap_score ( set1, set2 ):
    """ determines if one set overlaps with another ( 'semi-global jaccard' ) """
    shared = len( set1.__and__( set2 ) )
    return shared / float( min( len( set1 ), len( set2 ) ) )

def describe_median_rank ( m ):
    """ """
    cuts = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95, 1.0]
    names = ["very low", "low", "medium low", "trivial", "medium high", "high", "very high"]
    for cut, name in zip( cuts, names ):
        if m < cut:
            return name

# ---------------------------------------------------------------
# counting
# ---------------------------------------------------------------

counter = 0
kmer2term = {}
kmer2rank = {}
for items in reader( fh ):
    counter += 1
    term = items[0].strip()
    original = term
    if simplify:
        term = line.lower()
        term = re.sub( "[^a-z0-9]", " ", term )
        term = re.sub( " +", " ", term )
    kmers = term.split( " " ) if args.words else [term[j:j+k] for j in range( 0, len( term ) - k + 1 )]
    for kmer in kmers:
        temp = kmer2term.setdefault( kmer, {} )
        temp[original] = 1
        temp = kmer2rank.setdefault( kmer, {} )
        temp[counter] = 1

# ---------------------------------------------------------------
# identify the best kmers
# ---------------------------------------------------------------

medrank = {}
for kmer, ranksdict in kmer2rank.items():
    if len( ranksdict ) >= mincount:
        medrank[kmer] = 1 - median( ranksdict.keys() ) / float( counter )
limit = int( len( medrank ) * extent )
aSortedKmers = sorted( medrank.keys(), key=lambda x: medrank[x] )
aSortedKmers = aSortedKmers[0:limit] + aSortedKmers[-limit:]

# ---------------------------------------------------------------
# merge kmer clusters
# ---------------------------------------------------------------

aaKmerClusters = []
asetClusterTerms = []
for kmer in aSortedKmers:
    setTemp = set( kmer2term[kmer].keys() )
    # note: range( 0 ) returns [], so we execute the for..else either
    # ( i ) we failed to hit an existing cluster or ( ii ) the cluster list is empty ( base case )
    for i in range( len( aaKmerClusters ) ):
        if overlap_score( setTemp, asetClusterTerms[i] ) >= minmerge:
            aaKmerClusters[i].append( kmer )
            asetClusterTerms[i] = asetClusterTerms[i].__or__( setTemp )
            break
    else:
        aaKmerClusters.append( [kmer] )
        asetClusterTerms.append( setTemp )

# ---------------------------------------------------------------
# output to screen
# ---------------------------------------------------------------

print """Note:
MRANK if the average ( normalized ) rank of terms within the group
MRANK 1.0 is the exact top of the list
MRANK 0.5 is the exact middle
MRANK 0.0 is the exact bottom
"""

for i in range( len( aaKmerClusters ) ):
    # compute median median rank
    median_rank = median( [medrank[k] for k in aaKmerClusters[i]] )
    print "+----------------------------------------------"
    print "| GROUP : %05d" % ( i+1 )
    print "| MRANK :", median_rank, "( %s )" % ( describe_median_rank( median_rank ) )
    # print "KMERS :", "|".join( aaKmerClusters[i] )
    print "| ATOMS :"
    for kmer in aaKmerClusters[i]:
        print "\t", "<", kmer, ">"
    print "| TERMS :"
    for term in asetClusterTerms[i]:
        x = [k for k in aaKmerClusters[i] if k in term]
        #print "\t", term, "<--", "|".join( x )
        print "\t", term
    print
