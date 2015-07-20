#! /usr/bin/env python

import os, sys, re, glob, argparse
from zopy.table2 import table
from zopy.dictation import col2dict
from zopy.utils import path2name

dictMap = col2dict( sys.argv[1], key=0, value=1 )
tableData = table( sys.argv[2] )
tableData.apply_colheads( lambda x: path2name( x ) )
tableData.apply_colheads( lambda x: x.split( "." )[0] )
tableData.apply_colheads( lambda x: dictMap[x] )
tableData.dump( sys.argv[3] )
