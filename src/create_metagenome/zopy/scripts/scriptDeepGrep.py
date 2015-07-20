#! /usr/bin/env python

import os, sys, re, glob, argparse

"""
parser = argparse.ArgumentParser()
parser.add_argument( '-f', "--file_pattern" )
parser.add_argument( '-t', "--text_pattern" )
args = parser.parse_args()
file_pattern = args.file_pattern
text_pattern = args.test_pattern
"""

file_pattern = sys.argv[1]
text_pattern = sys.argv[2]

root = os.getcwd()
aMatchingFiles = []
aMatchingFiles += glob.glob( os.path.join( root, file_pattern ) )
for subroot, aDirectories, aFiles in os.walk( root ):
    for directory in aDirectories:
        aMatchingFiles += glob.glob( os.path.join( subroot, directory, file_pattern ) )

for name in aMatchingFiles:
    hit=False
    fh = open( name )
    for i, line in enumerate( fh ):
        if re.search( text_pattern, line ):
            if not hit:
                hit = True
                print name
            print "\t", "%5d" % ( i+1 ), line.strip()
    fh.close()
