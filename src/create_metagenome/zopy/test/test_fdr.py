#! /usr/bin/env python

import os, sys, re, glob, argparse
from random import random as rand
import zopy.fdr

alpha = 0.05
good_tests = 30
bad_tests = 70

p1 = [ ( rand() + 0.5 if i > good_tests/2 else rand() + 1.0 ) * alpha * i / ( good_tests + bad_tests ) for i in range( 1, good_tests+1 ) ]
p2 = [ rand() for i in range( bad_tests ) ]
pvalues = p1 + p2
qvalues = zopy.fdr.pvalues2qvalues( pvalues )
qvalues2 = zopy.fdr.pvalues2qvalues( pvalues, adjusted=True )
cutoff = zopy.fdr.fdr_cutoff( pvalues, alpha=alpha )
steps = zopy.fdr.steps( len( pvalues ), alpha=alpha )
test = [1 if p<=cutoff else 0 for p in pvalues]

# sort results
items = [( pvalues[i], qvalues[i], qvalues2[i], test[i] ) for i in range( len( pvalues ) )]
items = sorted( items )
for i in range( len( items ) ):
    print "\t".join( [str( k ) for k in list( items[i] ) + [steps[i]]] )
