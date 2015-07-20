#! /usr/bin/env python

import os, sys, re, glob, argparse
from zopy.enrichments import fisher_enrich as fe
from zopy.dictation import polymap
with open( sys.argv[1] ) as fh:
    x = [k.strip() for k in fh.readlines()]
a = polymap( sys.argv[2] )
for result in fe( x, a, ranked=True ):
    print result[:-1], len( result[-1] ), len( x )
