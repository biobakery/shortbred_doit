#! /usr/bin/env Rscript

#install.packages("optparse")
rm(list = ls())
library(optparse)
library(stringr)

#Function to read in the unix time files.
getData <- function(strName, strFile){
   sb1 <- file(strFile,"rt")
   strTime<- readLines(strFile,1)

   strTime <- str_match(strTime, "(.*)user")[,2]
   strResult = c(strTime)
   return(strResult)
}

option_list <- list(
  make_option("--sb1", default="sb1time.txt", type="character",
              help="Enter the first sb time file."),
  make_option("--cent1", default="centroid1time.txt", type="character",
              help="Enter the first cent time file."),
  make_option("--name1", default="ARDB", type="character",
              help="Enter the name of the first db."),
  make_option("--sb2", default="sb2time.txt", type="character",
              help="Enter the second sb time file."),
  make_option("--cent2", default="centroid2time.txt", type="character",
              help="Enter the second cent time file."),
  make_option("--name2", default="VF", type="character",
              help="Enter the name of the second db."),
  make_option("--reads", default=5000000, type="integer",
              help="Enter the number of reads for the synthetic metagenomes."),
  make_option("--out", default="Figure1.png", type="character",
              help="Enter the name of the output file.")
  
)
opt <- parse_args(OptionParser(option_list=option_list))


sb1<-opt$reads/as.numeric(getData(opt$name1,opt$sb1))
cent1<-opt$reads/as.numeric(getData(opt$name1,opt$cent1))
sb2<-opt$reads/as.numeric(getData(opt$name2,opt$sb2))
cent2<-opt$reads/as.numeric(getData(opt$name2,opt$cent2))

png(filename=opt$out,width=400)
barplot(height=c(sb1,cent1,sb2,cent2),col=c("red","blue"), 
        ylab="Reads per Second",
        names.arg=c(opt$name1,opt$name1,opt$name2,opt$name2))
legend("topright",c("ShortBRED","Centroids"),col=c("red","blue"),pch=15)
dev.off()


