#! /usr/bin/env Rscript
#source("http://bioconductor.org/biocLite.R")
#biocLite("genefilter")
rm(list=ls()) 
#install.packages("devtools")
#install.packages("multtest")
#install.packages("genomes")'multtest', 'genomes'
#library(devtools)
#install_github("LeviRmisc", user="lwaldron")
#writePCL() 

#Indicates which date to use in Jim's Dropbox for input/output files.
file.sep <- .Platform$file.sep

library(optparse)
#library(Biobase)
library(stringr)
#library(gplots)
#library(reshape)
#library(gridExtra)

dirARDB <-  paste("C:","Users","Jim","Dropbox","hutlab","Jim","ShortBRED","results","2013-09-16","metadata",sep=file.sep)
dirWD <- paste("C:","Users","Jim","Dropbox","hutlab","Jim","ShortBRED","results","2013-10-20","infant_metagenomes",sep=file.sep)

option_list <- list(
  make_option("--meta", default=paste(dirARDB,"ARDBmeta.txt",sep=file.sep), type="character",
              help="Enter the directory where the ARDB flat files stored."),
  make_option("--infant_results", default=paste(dirWD,"ARDB_InfantMetagenomesmergedresults.tab",sep=file.sep), type="character",
              help="Infant Results."),
  make_option("--by_gene", default=paste(dirWD,"ARDB_InfantMetagenomesmergedresults.tab",sep=file.sep), type="character",
              help="Infant Results."),
  make_option("--by_class", default=paste(dirWD,"ARDB_InfantMetagenomesmergedresults.tab",sep=file.sep), type="character",
              help="Infant Results.")
)
opt <- parse_args(OptionParser(option_list=option_list))


########################################################
# Load files
########################################################
dfARDB_Meta <- read.table(opt$meta, header=TRUE, sep="\t", row.names=1,
                  as.is=TRUE)
print(dim(dfARDB_Meta))
				  dfInfant_Results <- read.table(opt$infant_results, header=TRUE, sep="\t",row.names=1)
print(dim(dfInfant_Results))
dmInfant_Results<- data.matrix(dfInfant_Results)
print(dim(dmInfant_Results))
print("Loaded data...")


rownames(dmInfant_Results) <- rownames(dfInfant_Results)

#########################################################
print("Merging data...")
dfResults <- merge(dmInfant_Results,dfARDB_Meta,by="row.names")
print ("Data merged.")
print("Check first five column names...")
print(colnames(dfResults)[1:5])
print("Check last five column names...")
print(tail(colnames(dfResults), n=5))


print("Sum up by class...")
dfTest <- as.data.frame(aggregate(dfResults[,-which(names(dfResults) %in% c("Merged.ID","Row.names"))], by=list(dfResults$Merged.ID), FUN=sum))


filename <- opt$by_class
write.table(dfTest, file=filename, sep="\t")

filename <- opt$by_gene
write.table(dfResults, file=filename, sep="\t")
