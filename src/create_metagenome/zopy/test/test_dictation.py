#! /usr/bin/env python

import os, sys, re, glob, argparse
import zopy.dictation as d
file = sys.argv[1]

# ---------------------------------------------------------------
# prints
# ---------------------------------------------------------------

def print1 ( d ):
    for key, val in d.items():
        print key, "-->", val    
    print

def print2 ( d ):
    for key, d2 in d.items():
        print key
        for key2, val in d2.items():
            print "  ", key2, "-->", val
    print

# ---------------------------------------------------------------
# col2dict tests 
# ---------------------------------------------------------------

print "col2dict tests"
os.system( "head %s" % ( file ) )
print1( d.col2dict( file ) )
print1( d.col2dict( file, headers=True ) )
print1( d.col2dict( file, key=1, headers=True ) )
print1( d.col2dict( file, key=1, value=2, headers=True ) )
print1( d.col2dict( file, key=1, func=lambda row: float( row[2] ), headers=True ) )

# ---------------------------------------------------------------
# col2dict2 tests 
# ---------------------------------------------------------------

print "col2dict2 tests"
os.system( "head %s" % ( file ) )
print2( d.col2dict2( file ) )
print2( d.col2dict2( file, headers=True, mirror=True ) )
print2( d.col2dict2( file, headers=True, mirror=True, func=lambda row: float( row[2] ) ) )
print1( d.col2dict2( file, headers=True, func=lambda row: float( row[2] ), tupledict=True ) )
print1( d.col2dict2( file, headers=True, func=lambda row: float( row[2] ), tupledict=True, mirror=True ) )
