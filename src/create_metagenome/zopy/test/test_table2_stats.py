#! /usr/bin/env python

import os, sys, re, glob, argparse
from zopy.table2 import table
from zopy.table2_stats import one_sample, two_sample

t = table( sys.argv[1] )
t.head( 2, invert=True )
t.float()

aF1 = ["pielou"]
aF2 = ["bray_curtis", "spearman", "filtered_spearman"]

def sample ( rdict ):
    for key in sorted( rdict.keys(), key=lambda x: rdict[x], reverse=True ):
        print key, rdict[key]

for f in aF1:
    print "one sample", f
    sample( one_sample( t, f ) )

#for f in aF2:
for f in ["bray_curtis"]:
    print "two sample", f
    sample( two_sample( t, f ) )
