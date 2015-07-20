#!/usr/bin/env python
#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the SflE Scientific workFLow Environment for reproducible
# research, authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Curtis Huttenhower, chuttenh@hsph.harvard.edu).
#
# If you use this environment, the included scripts, or any related code in your work,
# please let us know, sign up for the SflE user's group (sfle-users@googlegroups.com),
# pass along any issues or feedback, and we'll let you know as soon as a formal citation
# is available.

import sys
import re

reTime = r'^(.*)user'
dictSB = {}
dictCent = {}

iCount = 1
for strLine in sys.stdin:
	if iCount%3==1:
		strMG, strType = strLine.split()
	elif iCount%3==2:
		mtchTime = re.search(reTime,strLine)
		dTime = float(mtchTime.group(1))
		if strType =="Centroid":
			dictCent[strMG]=dTime
		elif strType=="Marker":
			dictSB[strMG]=dTime
	iCount+=1

astrMGs = dictSB.keys()
astrMGs.sort()

print "\t".join(["Metagenome","SB Time","Centroid Time"])
for strMG in astrMGs:
	print "\t".join([strMG,str(dictSB[strMG]),str(dictCent[strMG])])

dMeanSB = sum(dictSB.values())/float(len(dictSB.values()))
dMeanCent = sum(dictCent.values())/float(len(dictCent.values()))

print "\t".join(["Mean:",'{0:.2f}'.format(dMeanSB),'{0:.2f}'.format(dMeanCent)])








