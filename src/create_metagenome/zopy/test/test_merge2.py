#! /usr/bin/env python

import os, sys, re, glob, argparse
from zopy.table2 import table

from time import time
start = time()
def interval ():
     global start
     stop = time()
     print stop - start, "seconds"
     start = stop

print "loading"
t1 = table( sys.argv[1] )
t2 = table( sys.argv[2] )
interval()
print "testing merge"
t1.merge( t2 )
interval()

print "loading"
t3 = table( sys.argv[1] )
t4 = table( sys.argv[2] )
interval()
print "testing merge_old"
t3.merge_old( t4 )
interval()

for i in range( len( t1.data ) ):
    for j in range( len( t1.data[0] ) ):
        if t1.data[i][j] != t3.data[i][j]:
            print "mismatch new:old", i, j, t1.data[i][j], t3.data[i][j]

for i in range( len( t3.data ) ):
    for j in range( len( t3.data[0] ) ):
        if t1.data[i][j] != t3.data[i][j]:
            print "mismatch old:new", i, j, t3.data[i][j], t1.data[i][j]
