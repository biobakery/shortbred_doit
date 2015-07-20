#! /usr/bin/env python

import os, sys, re, glob, argparse
from math import floor, log

dictLogCounts = { 0:0, "NA":0 }

for line in sys.stdin:
    try:
        n = float( line.strip() )
        if n == 0:
            dictLogCounts[0] += 1
        else:
            n = floor( log( n, 10 ) )
            n = "%.0e" % ( 10 ** n )
            if n not in dictLogCounts:
                dictLogCounts[n] = 0
            dictLogCounts[n] += 1
    except:
        dictLogCounts["NA"] += 1

total = sum( dictLogCounts.values() )
print " CUME   OMAG            :--------: 20%"
cume = 0
for logbin in sorted( dictLogCounts.keys(), key=lambda x: float( x ) if x != "NA" else -1 ):
    if logbin != "NA":
        value = dictLogCounts[logbin] / float( total )
        cume += value
        print "%.4f" % ( cume ), "%.4f" % ( value ), ":", logbin, "\t", "|" * int( 50 * value )
value = dictLogCounts["NA"] / float( total )
cume += value
print "%.4f" % ( cume ), "%.4f" % ( value ), ":", "NA", "\t", "|" * int( 50 * value )
print "total=", total
