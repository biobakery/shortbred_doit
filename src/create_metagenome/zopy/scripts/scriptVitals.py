#! /usr/bin/env python

import os, sys, re, glob, argparse
import csv
from numpy import mean, std
from scipy.stats.mstats import mquantiles

parser = argparse.ArgumentParser()
parser.add_argument( '-z', "--zeroes", action="store_true", help="if set: zeroes excluded" )
parser.add_argument( '-o', "--outliers", action="store_true", help="if set: outliers removed" )
args = parser.parse_args()

total = 0
bad = 0
aX = []
for aItems in csv.reader( sys.stdin, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE ):
    total += 1
    if len( aItems ) != 1:
       bad += 1
       print >>sys.stderr, "ignoring", "\t".join( aItems )
    else:
        aX.append( aItems[0] )

# determine sensible values
aX2 = []
for x in aX:
    try:
        x = float( x )
        if not args.zeroes or x != 0:
            aX2.append( x )
        else:
            bad += 1
    except:
        print >>sys.stderr, "ignoring", x
        bad += 1
aX = aX2

# output
q1, q2, q3 = mquantiles( aX )
iqr = q3 - q1


print "n   :", total
print "bad :", bad, "( %.1f%% )" % ( 100 * bad / float( total ) )
print "sum :", sum( aX )
print "min :", min( aX )
print "q1  :", q1
print "med :", q2
print "q3  :", q3
print "max :", max( aX )
print "avg :", mean( aX )
print "std :", std( aX )
print "cv  :", std( aX ) / mean( aX )
