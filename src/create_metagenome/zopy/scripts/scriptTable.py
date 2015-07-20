#! /usr/bin/env python

"""
Quick table manipulation
Uses unorthodox command line argument structure

Commands
--------
g: grep
h: head
n: normalize columns
s: select
t: transpose
u: keep only common features

Options
-------
c: by column
f: focus on non-header field (first arg)
i: invert
m: protect/reattach metadata
x: read choices from file (one per line)

===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os, sys
from zopy.table2 import table

args = sys.argv[1:]
# first arg is a command cluster; read stdin
if args[0][0] == "-":
    t = table()
# interpret first arg as a file name
else:
    t = table( args[0] )
    args = args[1:]

# determine operations
ops = []
i = 0
while i < len( args ):
    if args[i][0] == "-":
        op = [args[i], []]
        i += 1
        while i < len( args ) and args[i][0] != "-":
            op[1].append( args[i] )
            i += 1
        ops.append( op )

# hold metadata if asked
meta = None

# perform operations
for optype, opargs in ops:
    print >>sys.stderr,"COMMAND:", optype, opargs
    # process kwargs
    kwargs = {
        "transposed":True if "c" in optype else False,
        "invert": True if "i" in optype else False,
        }
    if "t" in optype:
        # transpose
        t.transpose()
    elif "m" in optype:
        # protect metadata
        if meta is None:
            meta, t = t.metasplit( opargs[0] )
        # restore metadata
        else:
            t.metamerge( meta )
            meta = None
    elif "n" in optype:
        # normalize columns
        t.float()
        t.normalize_columns()
        t.unfloat()
    elif "u" in optype:
        # unrarify
        t.float()
        t.unrarify( minlevel=float( opargs[0] ), mincount=float( opargs[1] ), **kwargs )
        t.unfloat()
    elif "h" in optype:
        # head (/tail, if inverted)
        t.head( opargs[0], **kwargs )
    # process args
    elif "s" in optype or "g" in optype:
        args = ["headers", opargs] if "f" not in optype else [opargs[0], opargs[1:]]
        if "x" in optype:
            with open( args[1][0] ) as fh:
                args[1] = [line.strip() for line in fh]
        t.select( *args, **kwargs ) if "s" in optype else t.grep( *args, **kwargs )

# reattach protected metadata
if meta is not None:
    t.metamerge( meta )

# return table
t.dump()
