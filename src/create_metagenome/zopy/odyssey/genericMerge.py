#! /usr/bin/env python

# ---------------------------------------------------------------
# convert Broad-style project output to single table
# ---------------------------------------------------------------

import os, sys, re, glob, argparse

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

default_fill_empty = "#N/A"
directory_sep = "/"

# ---------------------------------------------------------------
# argparse 
# ---------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument( "--input",          nargs="+", help="input file pattern (for glob)" )
parser.add_argument( "--key_col",        help="key column to slice from file (rowhead in final table)" )
parser.add_argument( "--value_col",      help="value column to slice from file" )
parser.add_argument( "--key_pattern",    help="pattern key (rowhead) must match, e.g. 'M[0-9]+:' for KEGG modules" )
parser.add_argument( "--fill_empty",     help="value to insert in place of missing values" )
parser.add_argument( "--strip_headers",  action="store_true", help="remove this file type's headers" )
parser.add_argument( "--strip_comments", action="store_true", help="strip comments from the file (lines beginning with #)" )
args = parser.parse_args()

# ---------------------------------------------------------------
# file loader 
# ---------------------------------------------------------------

def loader ( path, col1=0, col2=1, key_pattern=False, strip_headers=False, strip_comments=False ):
    fdict = {}
    fh = open(path)
    # headers?
    if strip_headers:
        fh.readline()
    # remaining lines
    for line in fh:
        items = line.strip().split("\t")
        # comment line?
        if strip_comments and items[0][0] == "#":
            continue
        # line too short?
        elif len(items) - 1 < max( col1, col2 ):
            continue
        # bad key?
        elif key_pattern and not re.search( key_pattern, items[col1] ):
            continue
        # include
        else:
            fdict[ items[col1] ] = items[col2]
    fh.close()
    return fdict

# ---------------------------------------------------------------
# discover all requested files; load their data
# ---------------------------------------------------------------

paths = args.input

# ---------------------------------------------------------------
# load all data 
# ---------------------------------------------------------------

data = {}
isfeature = {}
for path in paths:
    # get subfolder name
    fdict = loader( path, 
                    col1=int( args.key_col ), 
                    col2=int( args.value_col ), 
                    key_pattern=args.key_pattern, 
                    strip_headers=args.strip_headers, 
                    strip_comments=args.strip_comments 
                    )
    for f in fdict:
        isfeature[f] = 1
    if path not in data:
        data[ path ] = fdict
    else:
        print >>sys.stderr, "WARNING: TRYING TO ASSIGN 2ND FILE TO", path
        exit()

# ---------------------------------------------------------------
# make a table
# ---------------------------------------------------------------

# col and row headers
samples = sorted( data.keys() )
features = sorted( isfeature.keys() )

# set up col headers
outline = ["samples"] + [k for k in samples]
print "\t".join(outline)

# add rows
for feature in features:
    outline = [feature]
    for sample in samples:
        if feature in data[sample]:
            value = data[sample][feature]
        elif args.fill_empty:
            value = args.fill_empty
        else:
            value = default_fill_empty
        outline.append(value)
    print "\t".join(outline)

# ---------------------------------------------------------------
# report
# ---------------------------------------------------------------

print >>sys.stderr, "samples:", len(samples)
print >>sys.stderr, "features:", len(features)
