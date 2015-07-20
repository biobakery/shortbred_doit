#!/usr/bin/env python
"""
Created on Mon Oct  8 10:13:40 2012

@author: jim
"""

import sys
import argparse
import subprocess
import csv
import re
import random
import os

import Bio
from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers, and produces a relative abundance table.')


parser.add_argument('--nucs', type=str, dest='sNucs', help='Enter the path and name of the nucleotide file.')
parser.add_argument('--fastadir', type=str, dest='sFD', help='Enter the path to output the fasta files.', default ="fasta")
parser.add_argument('--genomes', type=str, dest='sGenomes', help='Enter the path and name of the file containing the names of the full genomes.')
parser.add_argument('--muG', type=float, dest='dGenomeMu', help='Mean of lognormal(mu,1) dist for genomes',default=1)
parser.add_argument('--muS', type=float, dest='dSpikeMu', help='Mean of lognormal(mu,1) dist for genomes',default=1)
parser.add_argument('-N', type=int, dest='iN', help='Enter the number of genes to include in the file.',default =10)
parser.add_argument('--padlength', type=int, dest='iPad', help='Enter the number of nucleotides to add to the gene at each end.',default =50)
parser.add_argument('--padgenome', type=str, dest='sPadGenome', help='Enter the path and name of the genome to use for padding.',default="pad.fasta")
parser.add_argument('--gold', type=str, dest='sGold', help='Enter the path and name of the gold standard file.',default="goldstandard.txt")
parser.add_argument('--abund', type=str, dest='sAF', help='Enter the path and name of the abundance file.',default="abundance.txt")
parser.add_argument('--dirgenomes', type=str, dest='sDirGenomes', help='Enter the path to the genomefiles.',default="")
parser.add_argument('--pctspike', type=float, dest='dPctSpike', help='Enter % reads that should come from spiked genes.',default=.05)
parser.add_argument('--log', type=str, dest='sLog', help='Enter the path and name of the log file.',default="log.txt")
parser.add_argument('--mutate', dest='bMutate', action='store_true', help='Add if you wish to mutate genes',default=False)

#parser.add_argument('--pct', type=float, dest='dPct', help='Enter the starting percentage, <= .30.', default=.20)
#parser.add_argument('--log', type=str, dest='sLog',help='Enter the name of the log file')
#parser.add_argument('--map', type=str, dest='sMap',help='Enter the name of the map file from usearch')


#parser.add_argument('--props', type=str, dest='sProps', help='Enter the path and name of the file with gene names and proportions.')

args = parser.parse_args()


dictNucs = {}
aaGoldStandard = []
iSetSize = 10
fileGS = open(args.sGold,'w')
aseqPad = []

def RandSlice(iLen, seqGene):
	#Pull a slice of iLen from seqGene, choose the starting nucleotide randomly
	iGeneLength = len(seqGene)
	iStart = random.randint(0,iGeneLength-(iLen+1))
	#print iGeneLength, iStart, (iStart+iLen)
	strPad = seqGene.seq[iStart:(iStart+iLen+1)]
	return strPad

def CountBases(strFileFasta):
	iBases = 0
	for seq in SeqIO.parse(strFileFasta, "fasta"):
		iBases = iBases + len(seq)
	return iBases

def check_create_dir( strDir ):

	if not os.path.exists( strDir ):
		os.makedirs( strDir )
	return strDir
 
def mutate_gene(seqIn, iMutations):
    dPoint = .60
    dInsert = .20
    #dDelete = 1 - dPoint - dInsert
    c_strNucs ="ATGC"
    
    for i in range(iMutations):
        dChoice = random.random()
        # Choose point mutation, insertion, or deletion
        if dChoice <= dPoint:
            iLoc = random.randint(0, len(seqIn)-1)
            strChange = c_strNucs[random.randint(0,len(seqIn)-1)]
            seqIn.seq[iLoc]= strChange
        elif dChoice > dPoint and dChoice <= (dPoint + dInsert):
            iLoc = random.randint(0, len(seqIn)-1)
            iLength = random.randint(1,6)
            strInsert = ""
            for i in range(iLength):
                strInsert = strInsert + c_strNucs[random.randint(0,len(seqIn)-1)]
            seqIn.seq = seqIn[0:iLoc] + strInsert + seqIn[iLoc:]
        elif dChoice > (dPoint + dInsert):
            iLength = random.randint(1,6)            
            iLoc = random.randint(0, len(seqIn)-1-iLength)
            seqIn.seq = seqIn[0:iLoc] + seqIn[iLoc:]
            
        
    

check_create_dir(args.sFD)

#Load in the padding material
for seq in SeqIO.parse(args.sPadGenome, "fasta"):
	aseqPad.append(seq)


#Load in the spiked genes, save them in dictNucs. dictNucs[SeqID] = sequence
for seq in SeqIO.parse(args.sNucs, "fasta"):
    seq.id = "USR_" + seq.id + "_END"
    dictNucs[seq.id] = seq

#Choose N spiked genes randomly, load in the names of the bugs/genomes from text file
setPickedGenes = random.sample(set(dictNucs.keys()), args.iN)
astrPickedGenomes = [line.strip() for line in open(args.sGenomes)]

#Break the N spiked genes into sets of (iSetSize), make a lognormal value for each set
iDraws = args.iN / iSetSize
if (args.iN % iSetSize >0):
	iDraws +=1

dRandArray = []
for i in range(iDraws):
    dRandArray = dRandArray + [(random.lognormvariate(args.dSpikeMu,1))]*iSetSize

#Give each of the bugs/genomes a lognormal value
#Make an array [BugName, LNDraw, Length, Length*LNDraw]
for strGeneName in astrPickedGenomes:
	dRand =random.lognormvariate(args.dGenomeMu,1)
	iBases = CountBases(args.sDirGenomes + os.sep + strGeneName)
	aGeneData = [strGeneName, dRand, iBases, iBases*dRand]
	aaGoldStandard.append(aGeneData)

#Sum up the total bug bases. we will set other parameters so that this is (1-dPctSpike)%
#of the bases of the metagenome. So if we want 5% of the reads to be spiked,
#we will make the bugs 95% of the metagenome (by BASES, not COUNT).

adBugBases = (zip(*aaGoldStandard)[3])
dTotBugBases = sum(adBugBases)

dPctBugs = float(1)-args.dPctSpike

dTotMGBases = dTotBugBases/dPctBugs


#Assign the lognormal draw for each spiked gene.
#Make an array where each row is [SpikedGeneName,lognormal val, GeneLength, TotGeneBases]
for strGeneName in setPickedGenes:
    aGeneData = [strGeneName, dRandArray.pop(), len(dictNucs[strGeneName]) + 2*args.iPad]
    aGeneData.append(aGeneData[1]*aGeneData[2])
    aaGoldStandard.append(aGeneData)

#Sum up the total bases for the spiked genes.
adSpikedBases = (zip(*aaGoldStandard[-(args.iN):])[3])
dTotSpikeBases = sum(adSpikedBases)

dFactor = (dTotMGBases*args.dPctSpike)/dTotSpikeBases

#Multiply each draw (for the spiked genes) by this factor
for aLine in aaGoldStandard[-(args.iN):]:
	aLine[1]=aLine[1]*dFactor
	aLine[3]=aLine[3]*dFactor



#Sum up the lognormal draws for everything (bugs/genomes + spiked genes), compute relative abundance
#Write this info to the goldstandard file.
adAbundance = (zip(*aaGoldStandard)[1])
dTotAbundance = sum(adAbundance)

adBases = (zip(*aaGoldStandard)[3])
dBases = sum(adBases)

#Print spiked genes as individual fasta files to "--fastadir"
aaGenes = aaGoldStandard[-(args.iN):]
for astrLine in aaGenes:
    dictNucs[astrLine[0]].seq = RandSlice(args.iPad,random.choice(aseqPad)) + dictNucs[astrLine[0]].seq + RandSlice(args.iPad,random.choice(aseqPad))
    SeqIO.write(dictNucs[astrLine[0]], args.sFD + os.sep + astrLine[0]+".fasta", "fasta")


fileGS.write("GeneName"+ "\t" + "Count" + "\t" + "Length" + "\t" + "TotBases" + "\t" +"RelCount" + "\t" + "RelBases" + "\n")
for aLine in aaGoldStandard:
    fileGS.write(aLine[0]+ "\t" + str(aLine[1]) + "\t" +str(aLine[2]) + "\t" + str(aLine[3]) +"\t" + str(aLine[1]/dTotAbundance) + "\t" + str(aLine[3]/dBases) + "\n" )





#Make the abundance text file
fileAbund = open(args.sAF,'w')


#Print lines for bugs/genomes
for aLine in aaGoldStandard[:-(args.iN)]:
	fileAbund.write(aLine[0] + "\t" + str("{0:.6f}".format(aLine[1]/dTotAbundance)) + "\n")


#Print lines for spiked genes
for aLine in aaGoldStandard[-(args.iN):]:
	fileAbund.write(aLine[0]+".fasta" + "\t" + str("{0:.6f}".format(aLine[1]/dTotAbundance)) + "\n")
	RandSlice(args.iPad,random.choice(aseqPad))

fileAbund.close()

"""
astrGeneRepeats = []

for aLine in aaGoldStandard:
    astrCopies =  aLine[1]*[aLine[0]]
    for x in astrCopies:
        astrGeneRepeats.append(x)

#print astrGeneRepeats
random.shuffle(astrGeneRepeats)
#print astrGeneRepeats


for strGene in astrGeneRepeats:
    SeqIO.write(dictNucs[strGene], sys.stdout, "fasta")

iTotBases = sum(zip(*aaGoldStandard)[3])
iTotGenes = sum(zip(*aaGoldStandard)[2])

fileGS.write("The shuffled and expanded genome file for " + args.sNucs + " contains" + "\n")

fileGS.write("Input genes: " + str(len(setPickedGenes)) + "\n")
fileGS.write("Gene copies: " + str(iTotGenes) + "\n")
fileGS.write("Total bases: " + str(iTotBases) + "\n")
"""
fileGS.close()
fileLog = open(args.sLog,'w')
fileLog.write("Ran fine.")
fileLog.close()



