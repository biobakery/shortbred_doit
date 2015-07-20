#!/usr/bin/env python

import os, sys, re, argparse

parser = argparse.ArgumentParser()
parser.add_argument( "template", help='' )
parser.add_argument( "-a", "--args", nargs="+", help='' )
parser.add_argument( "-e", "--execute", action="store_true", help='' )
args = parser.parse_args()

# utility
def path2name ( path ):
    items = os.path.split( path )[1].split( "." )
    return ".".join( items[0:len( items )-1] )

# load template
with open( args.template ) as fh:
    slurmlines = fh.readlines()

# make substitutions
for i, arg in enumerate( args.args ):
    i = i + 1
    slurmlines = map( 
        lambda line: re.sub( "\$%d([^0-9])" % ( i ), arg+"\\1", line ),
        slurmlines,
        )
    slurmlines = map( 
        lambda line: re.sub( "%%%d([^0-9])" % ( i ), path2name( arg )+"\\1", line ),
        slurmlines,
        )

# check that everything was substituted
for line in slurmlines:
    match = re.search( r"([$%][0-9]+)", line )
    if match:
        print >>sys.stderr, match.group( 1 ), "in", line, "not completed by arguments"

# write template
commandfile = "-".join( map( path2name, args.args ) ) + ".slurm"
with open( commandfile, "w" ) as fh:
    print >>fh, "".join( slurmlines )

# execute?
if args.execute:
    os.system( "sbatch %s" % ( commandfile ) )
    print >>sys.stderr, "submitted", commandfile
else:
    print >>sys.stderr, "Nothing executed; inspect slurm files for errors."
