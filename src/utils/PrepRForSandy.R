#! /usr/bin/env Rscript
rm(list = ls())
if(require("optparse")==FALSE){
	install.packages("optparse",lib="/data/jkaminski/Rpkgs",repos="http://cran.r-project.org")
	}
if(require("ROCR")==FALSE){
	install.packages("ROCR",lib="/data/jkaminski/Rpkgs",repos="http://cran.r-project.org")
	}
if(require("gridExtra")==FALSE){
  install.packages("gridExtra",lib="/data/jkaminski/Rpkgs",repos="http://cran.r-project.org")
}
if(require("reshape")==FALSE){
  install.packages("reshape",lib="/data/jkaminski/Rpkgs",repos="http://cran.r-project.org")
}
if(require("ggplot2")==FALSE){
  install.packages("ggplot2",lib="/data/jkaminski/Rpkgs",repos="http://cran.r-project.org")
}


library(optparse)
library(stringr)
library(ROCR)
