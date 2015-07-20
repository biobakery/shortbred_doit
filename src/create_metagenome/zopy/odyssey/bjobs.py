#! usr/bin/env python

import os, sys, re
counts = {}
for line in os.popen("bjobs"):
    line = line.strip()
    print >>sys.stderr, line
    items = re.split("\s+", line)
    if items[0] != "JOBID": # headers, weird lines
        if len(items) > 2:
            status = items[2]
            if status not in counts:
                counts[status] = 0
            counts[status] += 1
print "TOTAL JOBS =", sum(counts.values())
for status, count in counts.items():
    print "%s: %d" %(status, count),
print
