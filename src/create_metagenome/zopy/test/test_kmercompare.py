#! /usr/bin/env python

import os, sys, re, glob, argparse
import zopy.kmercompare as kmc

lines = [line.strip() for line in open( "../data/text.txt" )]
lines = [line for line in lines if line != ""]

n = len( lines )
for i in range( n ):
    for j in range( i+1, n ):
        a, b = kmc.compare( lines[i], lines[j] )
        #c, d = kmc.expectation( lines[i], lines[j] )
        if min( a, b ) > 0.7:# and max( c, d ) < 0.35:
            print i, j
            print lines[i]
            print lines[j]
            print a, b
            #print c, d
