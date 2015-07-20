#! usr/bin/env python

import random
cols = 1000
rows = 10000
outline = ["headers"]
for j in range( cols ):
    outline.append( "COL"+str( j+1 ) )
print "\t".join( outline )
for i in range( rows ):
    outline = ["ROW"+str( i+1 )]
    for j in range( cols ):
        outline.append( "%.5f" % ( random.random() ) if random.random() > 0.05 else "" )
    print "\t".join( outline )
