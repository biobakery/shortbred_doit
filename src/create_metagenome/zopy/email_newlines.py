#! usr/bin/env python

import re
import sys

def remove_newlines( filename ):
    fh = open( filename )
    string = ""
    for line in fh:
        string += line
    string = re.sub( "[\n\r]{2,}", "@@@@@", string )    
    string = re.sub( "[\n\r]", " ", string )    
    string = re.sub( " +", " ", string )
    string = re.sub( "@@@@@", "\n\n", string )
    print string

if __name__ == "__main__":
    remove_newlines( sys.argv[1] )
