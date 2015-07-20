#! usr/bin/env python

import random
cols = 2000
rows = 100000
outline = ["headers"]
for j in range( cols ):
    outline.append( "COL"+str( j+1 ) )
print "\t".join( outline )
"""
outline = ["META"]
for j in range( cols ):
    outline.append( random.choice( ["aaa", "bbb", "ccc"] ) )
print "\t".join( outline )
"""
for i in range( rows ):
    outline = ["ROW"+str( i+1 )]
    for j in range( cols ):
        #outline.append( "%.5f" % ( random.random() ) if random.random() > 0.05 else "" )
        outline.append( "%.5f" % ( random.random() ) )
        #outline.append( random.choice( ["ate", "bat", "dog", "man", "rig", "pig", "log", "fog", "imp", "dax", "roc", ""] ) )
    print "\t".join( outline )
