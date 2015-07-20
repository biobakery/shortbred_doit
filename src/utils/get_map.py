# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 10:50:56 2012

@author: jim
"""

import sys
import csv
import argparse

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers, and produces a relative abundance table.')


parser.add_argument('--uc', type=str, dest='sUC', help='Enter the path and name of the uc file')
parser.add_argument('--map', type=str, dest='sMap', help='Enter the path and name of the map file')
args = parser.parse_args()

def printMap(strUCMap,strTxtMap):
   
    dictGeneMap={}
    for strLine in csv.reader(open(strUCMap),delimiter='\t'):
        if (strLine[0] == "H"):
            dictGeneMap[strLine[-2]] = strLine[-1]
        elif (strLine[0] == "C"):
            dictGeneMap[strLine[-2]] = strLine[-2]           
            
    f = open(strTxtMap, 'w')
    for prot, fam in sorted(dictGeneMap.items(), key = lambda(prot, fam): (fam,prot)):
        f.write(fam + "\t" + prot + "\n")
    f.close()
    
printMap(args.sUC,args.sMap)
