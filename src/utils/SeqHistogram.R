#! /usr/bin/env Rscript

rm(list = ls())
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library("Biobase")

library("Biostrings")
library("optparse")

option_list <- list(
  make_option("--ff", default="C:\\Users\\Jim\\Desktop\\Histogram\\Markers20.faa", type="character",
              help="Enter a fasta file."),
  make_option("--wi", default=20, type="numeric",
              help="Enter the bin width."),
  make_option("--hst", default="C:\\Users\\Jim\\Desktop\\Histogram\\plot.png", type="character",
              help="Enter the output file.")
  
)

opt <- parse_args(OptionParser(option_list=option_list))

iBinWidth = opt$width
aseqDB <- readAAStringSet(opt$fasta, format="fasta",
                          nrec=-1L, skip=0L, use.names=TRUE)

aiLengths<-elementLengths(aseqDB)
bins=seq(0,max(aiLengths),by=iBinWidth)
png(opt$hplot)
hist(aiLengths,main="Histogram of Sequence Lengths",xlim=c(0,max(aiLengths)),xlab="",
     breaks=bins, sub=paste("Bin Width of",iBinWidth,",n=",length(aseqDB)))
dev.off()
