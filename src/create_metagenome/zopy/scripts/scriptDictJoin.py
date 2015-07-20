#! /usr/bin/env python

import os, sys, re, glob, argparse, csv
from collections import Counter

# argument parsing (python argparse)
parser = argparse.ArgumentParser()
parser.add_argument( 'file1', help='' )
parser.add_argument( 'file2', help='' )
parser.add_argument( '--focus1', default=1, type=int, help='base-1 column for phase1 key' )
parser.add_argument( '--focus2', default=1, type=int, help='base-1 column for phase2 key' )
parser.add_argument( '--remainder', action="store_true", help='print file2 lines that did not match' )
args = parser.parse_args()

c_sep = "\t"

# adjust to base-0 indexing for python
args.focus1 -= 1
args.focus2 -= 1

# load second file to dictionary
lengths2 = []
d = {}
with open( args.file2 ) as fh:
    for items in csv.reader( fh, dialect="excel-tab" ):
        lengths2.append( len( items ) )
        key = items[args.focus2]
        d.setdefault( key, {} )["\t".join( items )] = 1
print >>sys.stderr, "finished loading file2"

# make dummy line to add when join fails
if len( set( lengths2 ) ) != 1:
    exit( "file2 lines have unequal lengths" )
dummyline2 = "\t".join( "#N/A" for k in range( lengths2[0] ) )

# load second file, print join
counts = Counter()
lengths1 = []
hits = {}
with open( args.file1 ) as fh:
    for items in csv.reader( fh, dialect="excel-tab" ):
        line = "\t".join( items )
        lengths1.append( len( items ) )
        key = items[args.focus1]
        if key in d:
            hits[key] = 1
            counts[len( d[key] )] += 1
            for joinline in d[key]:
                print c_sep.join( [line, joinline] )
        else:
            counts[0] += 1
            print c_sep.join( [line, dummyline2] )

if args.remainder:
    if len( set( lengths1 ) ) != 1:
        sys.exit( "file1 lines have unequal lengths" )
    dummyline1 = "\t".join( ["#N/A" for k in range( lengths1[0] )] )
    for key in d:
        if key not in hits:
            for line in d[key]:
                print "\t".join( [dummyline1, line] )

# report 
print >>sys.stderr, "# summary"
for k in sorted( counts.keys() ):
    print >>sys.stderr, k, counts[k]
