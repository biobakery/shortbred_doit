#! /usr/bin/env python

"""
Align tabular data on the command line
-----------------------------------------------
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os, sys, re, glob, argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument( "-t", "--trim", action="store_true" )
parser.add_argument( "-n", "--norm", action="store_true" )
args = parser.parse_args()

c_iStringLimit = 30

aiMaxlens = []
aastrData = []
reader = csv.reader( sys.stdin, dialect="excel-tab" )
for astrRow in reader:
    if len( aiMaxlens ) == 0:
        aiMaxlens = [ 0 for strK in astrRow ]
    for iDex, strK in enumerate( astrRow ):
        if args.trim and len( strK ) >= c_iStringLimit:
            strK = strK[0:c_iStringLimit] + "..."           
            astrRow[iDex] = strK
        aiMaxlens[iDex] = max( aiMaxlens[iDex], len( strK ) )
    aastrData.append( astrRow )

if args.norm:
    iMaxMax = max( aiMaxlens )
    aiMaxlens = [iMaxMax] * len( aiMaxlens )

for astrRow in aastrData:
    for iDex, strK in enumerate( astrRow ):
        print strK + " " * ( aiMaxlens[iDex] - len( strK ) ), " ",
    print
