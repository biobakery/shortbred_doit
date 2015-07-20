#! /usr/bin/env python

"""
Slice columns from tables into dictionaries
"""

import sys
import csv

# ---------------------------------------------------------------
# shared methods
# ---------------------------------------------------------------

def funcGetValue ( row, value=None, func=None ):
    """ formats a value line for assignment to keys """
    if func and value:
        # get value, apply func, return
        return func( row[value] )
    elif value:
        # just return value
        return row[value]
    elif func:
        # apply function to whole row
        return func( row )
    else:
        # default value is true
        return True

def funcReadCSV ( path, headers=False ):
    """ shared file reader """
    fh = open( path )
    if headers:
        fh.readline()
    return csv.reader( fh, dialect="excel-tab" )

def funcTestInsert( dictD, key, value, verbose=True ):
    """ inserts a key, value pair but warns if we're overwriting something """
    if verbose and key in dictD and dictD[key] != value:
        print >>sys.stderr, "Warning: Overwriting", key, dictD[key], "with", value
    dictD[key] = value

def funcPair ( a, b ):
    """ Sort a pair of keys for smart tuplekey insertion """
    return ( a, b ) if a < b else ( b, a )

def funcMakeTupledict ( dictdictD, sort=True ):
    """
    Convert a nested dict [a][b] into a tupledict [( a, b )]
    default behavior is to perform standard sort on a, b to make keys unique
    """
    dictTemp = {}
    for key1, dictInner in dictdictD.items():
        for key2, value in dictInner.items():
            key = funcPair( key1, key2 ) if sort else ( key1, key2 )
            dictTemp[key] = value
    return dictTemp

# ---------------------------------------------------------------
# methods for loading dictionaries from files
# ---------------------------------------------------------------

def col2dict ( path, key=0, value=None, func=None, headers=False, verbose=True ):
    """
    From a tab-delimitted text file
    Return a nested dictionary where one columns is the KEY
    Additional VALUE column is optional
    """    
    dictD = {}
    for aRow in funcReadCSV( path, headers ):
        funcTestInsert( dictD, aRow[key], funcGetValue( aRow, value, func ), verbose=verbose )
    return dictD

def col2dict2 ( path, key1=0, key2=1, value=None, func=None, headers=False, mirror=False, tupledict=False, verbose=True ):
    """
    From a tab-delimitted text file
    Return a nested dictionary where two columns are the outer/inner key ( KEY1, KEY2 )
    Additional VALUE column is optional
    """
    dictdictD = {}
    for aRow in funcReadCSV( path, headers ):
        dictInner = dictdictD.setdefault( aRow[key1], {} )
        funcTestInsert( dictInner, aRow[key2], funcGetValue( aRow, value, func ), verbose=verbose )
        if mirror:
            dictInner = dictdictD.setdefault( aRow[key2], {} )
            funcTestInsert( dictInner, aRow[key1], funcGetValue( aRow, value, func ), verbose=verbose )
    return dictdictD if not tupledict else funcMakeTupledict( dictdictD, sort=mirror )

def polymap ( path, headers=False ):
    """
    input like:
    key1 val1 val2 val3
    key2 val4 val5
    output like:
    [key][value]
    """
    dictdictD = {}
    for aRow in funcReadCSV( path, headers ):
        dictInner = dictdictD.setdefault( aRow[0], {} )
        for value in aRow[1:]:
            dictInner[value] = True
    return dictdictD

# ---------------------------------------------------------------
# other methods
# ---------------------------------------------------------------

def join ( dictdictA, dictdictB ):
    """ [k1][k2] and [k2][k3]=value yields [k1][k3]=value """
    dictdictD = {}
    for k1, dictInnerA in dictdictA.items():
        for k2 in dictInnerA:
            for k3, value in dictdictB.setdefault( k2, {} ).items():
                dictdictD.setdefault( k1, {} )[k3] = value
    return dictdictD
