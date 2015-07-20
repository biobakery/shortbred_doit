#! /usr/bin/env python

# ---------------------------------------------------------------
# convert Broad-style project output to single table
# ---------------------------------------------------------------

import os, sys, re, glob, argparse
from zopy.table2 import table, c_strNA

# ---------------------------------------------------------------
# argparse 
# ---------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument( 
    "-1", "--input1", 
    help="first data table",
)
parser.add_argument( 
    "-2", "--input2", 
    help="second data table OR metadata table",
)
parser.add_argument( 
    "-o", "--output", 
    default=None, 
    help="output file",
)
parser.add_argument( 
    "-n", "--fill_empty", 
    default=None, 
    help="value to insert in place of missing values",
)
parser.add_argument( 
    "-m", "--metamerge", 
    action="store_true",
    help="specify if second table is a metadata table",
)
args = parser.parse_args()

# ---------------------------------------------------------------
# load and process data
# ---------------------------------------------------------------

table1 = table( args.input1 )
table2 = table( args.input2 )
if args.metamerge:
    table1.metamerge( table2 )
else:
    table1.merge( table2 )
if args.fill_empty is not None:
    table1.apply_entries( lambda x: x if x != c_strNA else args.fill_empty )

# ---------------------------------------------------------------
# dump table
# ---------------------------------------------------------------

# not ideal
if args.output is not None:
    table1.dump( args.output )
else:
    table1.dump()
