#! /usr/bin/env python

# ---------------------------------------------------------------
# convert Broad-style project output to single table
# ---------------------------------------------------------------

import os, sys, re, glob, argparse

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

directory_sep = "/"

# ---------------------------------------------------------------
# argparse 
# ---------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument( "--root",         help="root broad folder (contains per-sample folders)" )
parser.add_argument( "--file_pattern", help="pattern of desired files below root" )
parser.add_argument( "--col1",         help="key column to slice from file (rowhead in final table)" )
parser.add_argument( "--col2",         help="value column to slice from file" )
parser.add_argument( "--key_pattern",  help="pattern key (rowhead) must match, e.g. 'M[0-9]+:' for KEGG modules" )
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

paths = {}
for line in os.popen("find %s -name %s" % ( args.root, args.file_pattern )):
    paths[line.strip()] = 1

# ---------------------------------------------------------------
# load all data 
# ---------------------------------------------------------------

data = {}
isfeature = {}
for path in paths:
    # get subfolder name
    r = args.root
    r = r if r[-1] == directory_sep else r+directory_sep
    root_items = r.split( directory_sep )
    path_items = path.split( directory_sep )
    subfolder_name = path_items[ len(root_items)-1 ]
    fdict = loader( path, 
                    col1=int( args.col1 ), 
                    col2=int( args.col2 ), 
                    key_pattern=args.key_pattern, 
                    strip_headers=args.strip_headers, 
                    strip_comments=args.strip_comments 
                    )
    for f in fdict:
        isfeature[f] = 1
    if subfolder_name not in data:
        data[ subfolder_name ] = fdict
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
        value = "#N/A" if feature not in data[sample] else data[sample][feature]
        outline.append(value)
    print "\t".join(outline)

# ---------------------------------------------------------------
# report
# ---------------------------------------------------------------

print >>sys.stderr, "samples:", len(samples)
print >>sys.stderr, "features:", len(features)
