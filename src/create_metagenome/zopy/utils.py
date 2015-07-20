#! /usr/bin/env python

import os, sys, re, glob, argparse
import csv
from collections import defaultdict
from textwrap import fill

# ---------------------------------------------------------------
# text manipulation
# ---------------------------------------------------------------

def reader ( file_handle ):
    """ my favorite options for csv reader """
    for aItems in csv.reader( file_handle, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE ):
        yield aItems

def smartwrap ( text, charlim ):
    """ this function was already available in python... """
    return fill( text, width=charlim )

def path2name ( path, extgroups=1  ):
    """ given '/blah/blah/blah/hello.txt' returns 'hello'"""
    groups = os.path.split( path )[1].split( "." )
    return ".".join( groups[0:-extgroups] )

def tprint ( *args ):
    """ coerce list of items to strings then print with tabs between """
    print "\t".join( map( str, args ) )

def col2list ( filename, index=0, func=None, headers=False ):
    """ quickly load a column from a file """
    aItems = []
    with open( filename ) as fh:
        if headers:
            fh.readline()
        for aRow in reader( fh ):
            aItems.append( aRow[index] )
    return aItems if func is None else map( func, aItems )

def warn ( message ):
    print >>sys.stderr, "WARNING (%s):" % ( sys.argv[0] ), message

# ---------------------------------------------------------------
# dictionary methods
# ---------------------------------------------------------------

def sorteditems ( dictData, reverse=False, limit=None ):
    """ return k, v pairs in v-sorted order """
    iCounter = 0
    for key in sorted( dictData.keys(), key=lambda x: dictData[x], reverse=reverse ):
        iCounter += 1
        if limit is not None and iCounter >= limit:
            break
        else:
            yield ( key, dictData[key] )

def autodict ( iDepth=None, funcDefault=None ):
    """ 
    Acts as a constructor; makes an avdict 
    Can terminate at a specified depth as a defaultdict with the specified constructor.
    Example1: x = funcAVD( 3, int ) allows x["foo"]["bar"]["net"] += 1 ( i.e., a counter )
    Example2: x = funcAVD( 2, list ) allows x["foo"]["bar"].append( "net" ) 
    """
    if iDepth is None:
        return defaultdict( lambda: autodict( None ) )
    elif iDepth >= 2:
        return defaultdict( lambda: autodict( iDepth - 1, funcDefault ) )
    elif iDepth == 1:
        return defaultdict( funcDefault )

# ---------------------------------------------------------------
# tests
# ---------------------------------------------------------------

if __name__ == "__main__":
    import random
    x = autodict( 1, int )
    alpha = "abcdefghijklmnopqrstuvwxyz"
    for i in range( 1000 ):
        x[random.choice( alpha )] += 1
    for k, v in sorteditems( x, reverse=True, limit=10 ):
        print k, v
