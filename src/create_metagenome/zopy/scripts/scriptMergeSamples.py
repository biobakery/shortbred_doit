#! /usr/bin/env python

# ---------------------------------------------------------------
# convert Broad-style project output to single table
# ---------------------------------------------------------------

import os, sys, re, glob, argparse
from zopy.utils import reader, path2name
from zopy.table2 import nesteddict2table

# ---------------------------------------------------------------
# argparse 
# ---------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument( 
    "-i", "--input", 
    nargs="+", 
    help="input files",
)
parser.add_argument( 
    "-o", "--output",
    default=None,
    help="output file",
)
parser.add_argument( 
    "-k", "--key_col",
    default=0,
    help="key column",
)
parser.add_argument( 
    "-v", "--val_col",
    default=1,
    help="value column",
)
parser.add_argument( 
    "-p", "--key_pattern",
    default=None,
    help="pattern key ( rowhead ) must match, e.g. 'M[0-9]+:' for KEGG modules",
)
parser.add_argument( 
    "-n", "--fill_empty",
    default=None,
    help="value to insert in place of missing values",
)
parser.add_argument( 
    "-u", "--use_headers",  
    action="store_true", 
    help="use file headers as colheads and not file names",
)
parser.add_argument( 
    "-r", "--strip_headers",  
    action="store_true", 
    help="remove this file type's headers",
)
parser.add_argument( 
    "-c", "--strip_comments", 
    action="store_true", 
    help="strip comments from the file ( lines beginning with # )",
)
args = parser.parse_args()

# ---------------------------------------------------------------
# load all data
# ---------------------------------------------------------------

dictTableData = {}
# modified for faster looking up 4/2/2015
dictFeatureIndex = {}

for iDex, strPath in enumerate( args.input ):
	print >>sys.stderr, "Loading", iDex+1, "of", len( args.input ), ":", strPath
	aastrData = []
	strColhead = path2name( strPath )
	with open( strPath ) as fh:
		for astrItems in reader( fh ):
			aastrData.append( [astrItems[args.key_col], astrItems[args.val_col]] ),
	if args.strip_comments:
		aastrData = [ astrRow for astrRow in aastrData if astrRow[0][0] != "#" ]
	if args.use_headers:
		strColhead = aastrData[0][1]
	if args.strip_headers:
		aastrData = aastrData[1:]
	if args.key_pattern:
		aastrData = [ astrRow for astrRow in aastrData if re.search( args.key_pattern, astrRow[0] ) ]
	for strFeature, strValue in aastrData:
		if strFeature not in dictFeatureIndex:
			dictFeatureIndex[strFeature] = len( dictFeatureIndex ) + 1
	dictTableData[strColhead] = { strKey:strValue for [strKey, strValue] in aastrData }

# ---------------------------------------------------------------
# coerce to table
# ---------------------------------------------------------------

# not ideal
kwargs={"empty":args.fill_empty} if args.fill_empty is not None else {}

# try to maintain original ordering (modified 4/2/2015)
astrFeatures = sorted( dictFeatureIndex.keys(), key=lambda x: dictFeatureIndex[x] )

tableData = nesteddict2table( dictTableData, aColheads=astrFeatures, **kwargs )
tableData.rowsort()
tableData.transpose()

# not ideal
if args.output is not None:
    tableData.dump( args.output )
else:
    tableData.dump()
