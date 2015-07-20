#! usr/bin/env python

import sys
from zopy.table2 import table
t = table( sys.argv[1] )
t.float()
t.groupby( lambda x: x[0].split( "|" )[0], sum )
t.dump()
