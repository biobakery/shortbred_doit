#! /usr/bin/env python

import sys
from Bio import SeqIO
minLength = None
for record in SeqIO.parse( sys.stdin, "fastq" ):
    myLength = len( record.seq )
    if not minLength or myLength < minLength:
        minLength = myLength
        print minLength

