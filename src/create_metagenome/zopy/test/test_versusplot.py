#! usr/bin/env python

from versusplot import versusplot
import random

# generate two vectors of roughly correlated random data
w = [1 + 5 * random.random() for i in range( 500 )]
w = [10**-k for k in w]
x = [0.0 if random.random() < 0.25 else 1e-6 + k * random.random() for k in w]
y = [0.0 if random.random() < 0.75 else 1e-6 + k * random.random() for k in w]
# isolate a subset of "special" points
col = []
zorder = []
pch = []
groups = []
choices = [( "1", "orange" ), ( "2", "green" ), ( "3", "purple" )]
for k in x:
    if random.random() < 0.1:
        choice = random.choice( choices )
        col.append( choice[1] )
        groups.append( choice[0] )
        zorder.append( 1 )
        pch.append( 4 )
    else:
        col.append( "#00000025" )
        groups.append( "0" )
        zorder.append( 0 )
        pch.append( 16 )
versusplot( x, y, xlabel="random x", ylabel="random y", title="my plot", logmin=-6, logmax=-1, col=col, zorder=zorder, pch=pch, groups=groups )
