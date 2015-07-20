#! /usr/bin/env Rscript

#install.packages("optparse")
rm(list = ls())
library(optparse)

option_list <- list(
#  make_option(c("-i", "--inputdir"), default="C:\\Users\\Jim\\Dropbox\\hutlab\\Jim\\thesis\\results\\eval2010-12-05\\",
#              help="Enter the directory containing the csv files.")
  make_option(c("-N", "--number"), default=10,
              help="Enter a number to be multiplied by 3")  
)
opt <- parse_args(OptionParser(option_list=option_list))
print(opt$number)
