#! usr/bin/env python

#import table
from table2 import table

print "method 1"
m = table( "big1.txt" )
f = table( "big2.txt" )
m.augment( f )
print "method 2"
m = table( "big1.txt" )
f = table( "big2.txt" )
m.augment2( f )



"""
data = [["%d.%d" % ( x,y ) for x in range( 10 )] for y in range( 20 )]
t = table2.table( data )
t.dump( "temp2.txt" )
"""

"""
t.grep( 0, "ROW1" )
t2 = t.grep( 0, "ROW.*0$", in_place=False )
t.dump( "temp2.txt" )
t2.dump( "temp3.txt" )
print t.size()
print t2.size()
"""

"""
t.stratify( "META" )
"""

"""
metadata = table2.table( "temp.txt" )
metadata.select( "headers", "META" )
metadata.delete( "headers", "COL5", transposed=True )

data = table2.table( "temp.txt" )
data.grep( "headers", "2|4|6|8", transposed=True )

metadata.augment( data )
metadata.dump( "temp2.txt" )
"""
