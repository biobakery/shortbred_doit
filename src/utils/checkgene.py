#!/usr/bin/env python
import argparse
import re
import csv
import sys
import string

import Bio
#from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Check evaluation results.')

#Input
parser.add_argument('--markers', type=str, dest='sMarkers', help='Enter the path and name of the SB markers faa file).')
parser.add_argument('--centroids', type=str, dest='sCentroids', help='Enter the path and name of the centroids faa file.')
parser.add_argument('--wgs', type=str, dest='sWGS', help='Enter the path and name of the wgs.fna.')
parser.add_argument('--markerblast', type=str, dest='sMarkerBlast', help='Enter the path and name of the blast results for the SB markers.')
parser.add_argument('--centblast', type=str, dest='sCentBlast', help='Enter the path and name of the blast results for the centroids')
parser.add_argument('--qmarkers', type=str, dest='sMarkerSB', help='Enter the path and name of the SBQuantify results for the markers.')
parser.add_argument('--qcent', type=str, dest='sCentSB', help='Enter the path and name of the SBQuantify results for the centroids')
parser.add_argument('--map', type=str, dest='sMap', help='Enter name and path of the map file')
parser.add_argument('--debug', type=str, dest='sDebug', default="",help='Enter a gene to check in detail')
parser.add_argument('--nucs', type=str, dest='sNucs', default="",help='Enter the file containing nucleotides')
parser.add_argument('--FPfile', type=str, dest='sFPfile', default="",help='Enter where the FPs will print')
parser.add_argument('--reads', type=int, dest='iReads', default=5000000,help='Enter the number of reads in the synthetic metagenome')


#Parameters
parser.add_argument('--tmid', type=float, dest='dTMID', help='Enter the percent identity for a TM match', default = .95)
parser.add_argument('--qmid', type=float, dest='dQMID', help='Enter the percent identity for a QM match', default = .95)
parser.add_argument('--alnlength', type=int, dest='iAlnLength', help='Enter the minimum alignment length. The default is 20', default = 20)
parser.add_argument('--alnTM', type=int, dest='iAlnMax', help='Enter a bound for TM alignments, such that aln must be>= min(markerlength,alnTM)', default = 25)

args = parser.parse_args()

astrGenes =[]
dictFams={}

dictCount = {}
dictSBMarkers = {}
dictSBMarkerLength={}
dictSBCent = {}

dictLengths = {}

iReads = args.iReads

f = open(args.sFPfile,'w')

for astrLine in csv.reader( open(args.sMap), csv.excel_tab ):
            dictFams[str(astrLine[1]).strip()]=str(astrLine[0]).strip()

for geneNuc in SeqIO.parse(args.sNucs, "fasta"):
    dictLengths[geneNuc.id] = len(geneNuc.seq)



#Load WGS count
#This a tab-delimited file with the form (gene,count,family)
for astrLine in csv.reader( open(args.sWGS, 'rb'), delimiter='\t' ):
    strFam= astrLine[2].strip()
    dictCount[strFam] = dictCount.get(strFam,0) + int(astrLine[1])

#Load Quantify output for centroids
#Adjusted on 3/25/2013, new ShortBRED-Quantify output has a one line header.
#The old version had 12 lines of usearch stderr output at the beginning.
iLineCount=0
for astrLine in csv.reader( open(args.sMarkerSB, 'rb'), delimiter='\t' ):
    #print astrLine
    iLineCount+=1
    if (iLineCount>1):
        dictSBMarkers[astrLine[0]] = astrLine[1]
        dictSBMarkerLength[astrLine[0]] = astrLine[3]

#Load Quantify output for markers
#Adjusted on 3/25/2013, new ShortBRED-Quantify output has a one line header.
#The old version had 12 lines of usearch stderr output at the beginning.
iLineCount=0
for astrLine in csv.reader( open(args.sCentSB, 'rb'), delimiter='\t' ):
    iLineCount+=1
    if (iLineCount>1):
        dictSBCent[astrLine[0]] = astrLine[1]


#Get Marker Count
dictTMs = {}
dictQMs = {}

for strLine in open(args.sMarkers, 'rb'):
    mtchTM = re.search(r'>(.*)_TM.*',strLine)
    if(mtchTM):
        strStub = mtchTM.group(1)
        dictTMs[strStub] = dictTMs.get(strStub,0) +1

	# I will treat the JM's the same as QM's for this part.
    mtchQM = re.search(r'>(.*)_[QJ]M.*',strLine)
    if(mtchQM):
        strStub = mtchQM.group(1)
        dictQMs[strStub] = dictQMs.get(strStub,0) +1


#Get Marker Lengths
dictMarkerLen = {}
dictMarkerLenAll = {}


for seq in SeqIO.parse(args.sMarkers, "fasta"):
 	dictMarkerLen[seq.id] = len(seq)

#Get Count in ShortBRED Blast Results

dictMarkerHits = {}
dictMarkerTP = {}
dictMarkerFP = {}

iLine=0
for astrLine in csv.reader( open(args.sMarkerBlast, 'rb'), delimiter='\t' ):
    #print(astrLine[1])
    mtchTM = re.search(r'_TM.*',astrLine[1])

    if (mtchTM):
        dID = args.dTMID
        iAln = min(dictMarkerLen[astrLine[1]] ,args.iAlnMax)
    else:
        dID = args.dQMID
        iAln = args.iAlnLength

    if (int(astrLine[3])>= iAln and (float(astrLine[2])/100.0) >= dID):
        mtchGene = re.search(r'USR_(.*)_END',astrLine[0])
        mtchMarkerFam = re.search(r'(.*)_[TQJ]M.*',astrLine[1])
        strMarkerFam = mtchMarkerFam.group(1)
        dictMarkerHits[strMarkerFam] = dictMarkerHits.get(strMarkerFam,0) +1

        if(mtchGene):
            strGeneHit = str(mtchGene.group(1)).strip()
            if (dictFams.get(strGeneHit,"NotInDict")==dictFams[strMarkerFam]):
                dictMarkerTP[strMarkerFam] = dictMarkerTP.get(strMarkerFam,0) +1
            else:
                dictMarkerFP[strMarkerFam] = dictMarkerFP.get(strMarkerFam,0) +1
                f.write("Marker FP:" + ('\t'.join(astrLine))+ "\n")
        else:
            dictMarkerFP[strMarkerFam] = dictMarkerFP.get(strMarkerFam,0) +1
            f.write("Marker FP:" + ('\t'.join(astrLine))+ "\n")


#get count in Centroid results
dictCentHits = {}
dictCentTP = {}
dictCentFP = {}

for astrLine in csv.reader( open(args.sCentBlast, 'rb'), delimiter='\t' ):
    mtchGene = re.search(r'USR_(.*)_END',astrLine[0])
    strCentFam = astrLine[1]
    dictCentHits[strCentFam] = dictCentHits.get(strCentFam,0) +1
    if(mtchGene):
        strGeneHit = str(mtchGene.group(1)).strip()
        if (dictFams.get(strGeneHit,"NotInDict")==dictFams[strCentFam]):
            dictCentTP[strCentFam] = dictCentTP.get(strCentFam,0) +1
        else:
            dictCentFP[strCentFam] = dictCentFP.get(strCentFam,0) +1
    else:
        dictCentFP[strCentFam] = dictCentFP.get(strCentFam,0) +1
"""
iMarkerFN = int(dictCount.get(strName,0)) - iMarkerTP
iMarkerTN = iReads - iMarkerTP - iMarkerFP - iMarkerFN
"""

print "\t".join(["GeneName","WGS Reads","TMs","QMs","Q-Markers","Q-Centroids","MarkerHits","MarkerTP","MarkerFP","MarkerTN","MarkerFN","CentroidHits","CentroidTP","CentroidFP","CentroidTN","CentroidFN","TotMarkerLength","NucLength"])

setGenes = set(dictFams.values())


for strName in setGenes:
    #Get WGS Count
    iWGSCount =  str( dictCount.get(strName,0))
    iTM = dictTMs.get(strName,0)
    iQM = dictQMs.get(strName,0)
    iMarkerHits = dictMarkerHits.get(strName,0)
    iMarkerTP = dictMarkerTP.get(strName,0)
    iMarkerFP =dictMarkerFP.get(strName,0)
    iMarkerFN = int(dictCount.get(strName,0)) - iMarkerTP
    iMarkerTN = iReads - iMarkerTP - iMarkerFP - iMarkerFN

    iCentHits =dictCentHits.get(strName,0)
    iCentTP = dictCentTP.get(strName,0)
    iCentFP = dictCentFP.get(strName,0)
    iCentFN = int(dictCount.get(strName,0)) - iCentTP
    iCentTN = iReads - iCentTP - iCentFP - iCentFN
    iTotMarkerLen = dictSBMarkerLength.get(strName,0)
    iNucLen = dictLengths.get(strName,0)


    print "\t".join([strName,str(iWGSCount),str(iTM),str(iQM),str(dictSBMarkers.get(strName,0)),
	str(dictSBCent.get(strName,0)),str(iMarkerHits),str(iMarkerTP),str(iMarkerFP),str(iMarkerTN),
	str(iMarkerFN),str(iCentHits),str(iCentTP),str(iCentFP),str(iCentTN),str(iCentFN),str(iTotMarkerLen),str(iNucLen)])


"""
if (args.sDebug!=""):
    strName = str(args.sDebug)

    iWGSCount =  str( dictCount.get(strName,0))


    #Get Marker Count

    iTM=0
    iQM=0
    for strLine in open(args.sMarkers, 'rb'):
        mtchTM = re.search(r'>'+ re.escape(strName) + r'_TM_.*',strLine)
        if(mtchTM):
            iTM+=1
            print "TM: " + strLine
        mtchQM = re.search(r'>'+ re.escape(strName) + r'_QM_.*',strLine)
        if(mtchQM):
            iQM+=1
            print "QM: " + strLine



    #Get Count in ShortBRED Blast Results
    iMarkerHits = 0
    iMarkerTP = 0
    iMarkerFP = 0
    for astrLine in csv.reader( open(args.sMarkerBlast, 'rb'), delimiter='\t' ):
        mtch = re.search(re.escape(strName),astrLine[1])
        if mtch:
            iMarkerHits += 1
            mtchGene = re.search(r'USR_(.*)_END',astrLine[0])
            if(mtchGene):
                strGeneHit = str(mtchGene.group(1)).strip()
                if (dictFams[strGeneHit]==dictFams[strName]):
                    iMarkerTP+=1
                    print "Marker TP:" + str(astrLine[0]) +"\t" + str(astrLine[1])
                else:
                    iMarkerFP+=1
                    print "Marker FP:" + str(astrLine[0]) +"\t" + str(astrLine[1])
            else:
                iMarkerFP+=1
                print "Marker FP:" + str(astrLine[0]) +"\t" + str(astrLine[1])
    iMarkerFN = int(dictCount.get(strName,0)) - iMarkerTP
    iMarkerTN = iReads - iMarkerTP - iMarkerFP - iMarkerFN



    #print strName, iMarkerHits


    #Get Count in Centroid Blast Results
    iCentHits = 0
    iCentTP = 0
    iCentFP = 0
    for astrLine in csv.reader( open(args.sCentBlast, 'rb'), delimiter='\t' ):
        mtch = re.search(re.escape(strName),astrLine[1])
        if mtch:
            iCentHits += 1
            mtchGene = re.search(r'USR_(.*)_END',astrLine[0])
            if(mtchGene):
                strGeneHit = str(mtchGene.group(1)).strip()
                if (dictFams[strGeneHit]==dictFams[strName]):
                    iCentTP+=1
                    print "Cent TP:" + str(astrLine[0]) +"\t" + str(astrLine[1])
                else:
                    iCentFP+=1
                    print "Cent FP:" + str(astrLine[0]) +"\t" + str(astrLine[1])
            else:
                iCentFP+=1
                print "Cent FP:" + str(astrLine[0]) +"\t" + str(astrLine[1])
    iCentFN = int(dictCount.get(strName,0)) - iCentTP
    iCentTN = iReads - iCentTP - iCentFP - iCentFN
"""
