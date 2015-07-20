#! /usr/bin/env python

import os, sys, re, glob, argparse
import csv
from scipy.stats import spearmanr
from scipy.stats.mstats import mquantiles
from zopy.stats import mutinfo

parser = argparse.ArgumentParser()
parser.add_argument( '-z', "--zerozero", action="store_true", help="if set: 0,0 points EXCLUDED" )
parser.add_argument( '-q', "--quantmi", action="store_true", help="if set: quantize and do normalized mutual information" )
args = parser.parse_args()

total = 0
bad = 0
aX = []
aY = []
for aItems in csv.reader( sys.stdin, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE ):
    total += 1
    if len( aItems ) != 2:
       bad += 1
       print >>sys.stderr, "ignoring", "\t".join( aItems )
    else:
        aX.append( aItems[0] )
        aY.append( aItems[1] )

# determine sensible pairs
aX2 = []
aY2 = []
for x, y in zip( aX, aY ):
    try:
        x = float( x )
        y = float( y )
        if max( x, y ) > 0 or not args.zerozero:
            aX2.append( x )
            aY2.append( y )
        else:
            bad += 1
    except:
        print >>sys.stderr, "ignoring", x, y
        bad += 1
aX, aY = aX2, aY2

# quantile mutinfo
def quantform ( numbers ):
    quantiles = mquantiles( numbers )
    output = []
    for n in numbers:
        for i, q in enumerate( quantiles ):
            if n <= q:
                output.append( "Q%d" % ( i+1 ) )
                break
        else:
            output.append( "Q%d" % ( 1 + len( quantiles ) ) )
    return output

# output
r, p = spearmanr( aX2, aY2 )
print "total pairs  :", total
print "bad pairs    :", bad, "( %.1f%% )" % ( 100 * bad / float( total ) )
print "spearman     :", r
print "p-value      :", p
print "norm mutinfo :", "-" if not args.quantmi else mutinfo( quantform( aX ), quantform( aY ), normalized=True )
