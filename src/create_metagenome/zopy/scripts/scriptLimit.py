#! /usr/bin/env python

import os, sys, re, glob, argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument( '-c', "--col",   type=int,     help="focal column, base-1 indexing" )
parser.add_argument( '-l', "--lower", default=None, type=float, help="lower limit ( none if undef )" )
parser.add_argument( '-u', "--upper", default=None, type=float, help="upper limit ( none if undef )" )
args = parser.parse_args()

for aItems in csv.reader( sys.stdin, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE ):
    try:
        value = float( aItems[args.col - 1] )
        include = True
        if args.lower is not None and ( value < args.lower ):
            include = False
        if args.upper is not None and ( value > args.upper ):
            include = False
        if include:
            print "\t".join( aItems )
    except:
        continue
