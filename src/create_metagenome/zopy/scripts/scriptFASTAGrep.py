#! /usr/bin/env python

import sys, re, argparse
from Bio import SeqIO

# argument parsing (python argparse)
parser = argparse.ArgumentParser()
parser.add_argument( 
    '-i', '--input', 
    type=argparse.FileType('r'), 
    default=sys.stdin, 
    help='fasta file', )
parser.add_argument( 
    '-p', '--pattern', 
    default="", 
    help='perl-like regex pattern', )
parser.add_argument( 
    '-v', '--inverted', 
    action="store_true", 
    help='inverted search', )
parser.add_argument( 
    '-t', '--tabbed', 
    action="store_true", 
    help='tabbed output', )
args = parser.parse_args()

# search
for record in SeqIO.parse( args.input, "fasta" ) :
    match = re.search( args.pattern, record.name )
    if ( match and not args.inverted ) or ( not match and args.inverted ):
        if args.tabbed:
            print "\t".join( [record.name, str( record.seq )] )
        else:
            print ">"+record.name
            print str( record.seq )
