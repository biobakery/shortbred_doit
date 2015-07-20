#! /usr/bin/env Rscript

#install.packages("optparse",lib="/home/jkaminski/Rpkgs",repos="http://cran.r-project.org")
rm(list = ls())
.libPaths('/n/home11/jkaminski/Rlibs')
library(optparse)
library(stringr)

option_list <- list(
  make_option("--gs", default="gs.txt", type="character",
              help="Enter the gold standard file."),
  make_option("--wc", default="wgscounts.txt", type="character",
              help="Enter the file of wgs counts."),
  make_option("--out", default="output.txt", type="character",
              help="Enter the name of the output file."),
  make_option("--main", type="character", default="",
              help="Enter the name of the main output file."),
  make_option("--high", type="character", default="",
              help="Enter the output file for high families."),  
  make_option("--low", type="character", default="",
              help="Enter the output file for low families.") , 
  make_option("--cent", type="character", default="",
              help="Enter the length of the centroids."),
  make_option("--markers", type="character", default="",
              help="Enter marker results file.")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt$gs)


#Load GoldStandard 
dfGS <- read.table(opt$gs,header=TRUE,sep="\t")
dfGS[,1] <- str_match(dfGS[,1], "USR_(.*)_END")[,2]
write.table(dfGS, "outGS.txt", sep="\t", quote=FALSE)

#Load WGS Counts
dfWGS <- read.table(opt$wc,header=FALSE,sep="\t")
colnames(dfWGS)<- c("GeneName","Reads","Family")

dfMerged <- merge(dfGS,dfWGS,by="GeneName",all=TRUE,sort=FALSE)
write.table(dfMerged, "outMerged.txt", sep="\t", quote=FALSE)


#Sum up by family - Following Curtis's suggestion, I am now using "RelCount", 
# which is GeneCount/TotalGenes in Metagenome. Previously used "Count". 4/19/2013

dfCountByFamily <- aggregate.data.frame(cbind(dfMerged$RelCount,dfMerged$RelBases),by=list(dfMerged$Family),sum)
colnames(dfCountByFamily)<- c("Family","RelCount","RelBases")
write.table(dfCountByFamily, "outFam.txt", sep="\t", quote=FALSE)


#Merge Family Counts onto the main file.
dfMain <- read.table(opt$main,header=TRUE,sep="\t")
colnames(dfMain)[1] <- "Family"
dfMainFams <- merge(dfMain,dfCountByFamily,by="Family",all=TRUE,sort=FALSE)
#write.table(dfGS, "outGS.txt", sep="\t", quote=FALSE)

#Load the centroid lengths
dfCentroids <- read.table(opt$cent,header=FALSE,sep=" ")
colnames(dfCentroids)<- c("Family","CentLength")
dfMainFams <- merge(dfMainFams,dfCentroids,by="Family",all.x=TRUE,sort=FALSE)

dfMainFams$SBRelCount <- dfMainFams$Q.Markers/sum(dfMainFams$Q.Markers)

dfMainFams$SBRelCount[is.na(dfMainFams$SBRelCount)] <- 0
dfMainFams$RelCount[is.na(dfMainFams$RelCount)] <-0
dfMainFams$RelBases[is.na(dfMainFams$RelBases)] <-0

dfMainFams$High <- 0
dfMainFams$High[ (( (dfMainFams$SBRelCount / dfMainFams$RelCount) > 1.2) | (dfMainFams$SBRelCount > 0 & dfMainFams$RelCount == 0))] <- 1

dfMainFams$Low <- 0
dfMainFams$Low[ (( (dfMainFams$SBRelCount / dfMainFams$RelCount) < 0.8) | (dfMainFams$SBRelCount == 0 & dfMainFams$RelCount > 0))] <- 1

write.table(dfMainFams, opt$out, sep="\t", quote=FALSE,row.names=FALSE)

print("Made eval table")
#Make a table of basic information for the Low/High Families.
dfMainFams <- dfMainFams[, c(1:5,7:9,19,22,23,24) ]

#Merge on the marker level data
dfMarkers <- read.delim(opt$markers,header=TRUE,sep="\t")
print("Loaded markers")
dfMainFams <- merge(dfMarkers,dfMainFams,by="Family",all.y=TRUE,sort=FALSE)

dfHigh <- dfMainFams[dfMainFams$High==1,]
dfLow <- dfMainFams[dfMainFams$Low==1,]

write.table(dfHigh, opt$high, sep="\t", quote=FALSE,row.names=FALSE)
write.table(dfLow, opt$low, sep="\t", quote=FALSE,row.names=FALSE)
