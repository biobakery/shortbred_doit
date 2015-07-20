#! /usr/bin/env Rscript

#install.packages("optparse")
rm(list = ls())
library(optparse)
library(stringr)
library(ROCR)

option_list <- list(
  make_option("--AUC", default="C:\\Users\\Jim\\Desktop\\VFallAUC.txt", type="character",
              help="Enter the name of the file containing AUC values"),
  make_option("--db", default="VF", type="character",
              help="Enter the db used."),
  make_option("--outdir", default="", type="character",
              help="Enter the directory for output.")
)

opt <- parse_args(OptionParser(option_list=option_list))
dfResults <- read.table(opt$AUC,header=TRUE,sep="\t")
#colnames(dfResults) <- c("DB","PctID","MarkerLen","MatchID","SB_AUC","C_AUC","SB_SCorr","C_SCorr","SB_Spec","C_Spec","SB_Sens","C_Sens")

dfResults <- dfResults[dfResults$Dataset==opt$db,]
print(dim(dfResults))
#attach(dfResults)

#Plot the AUC for ShortBRED

listMatchIDs = levels(as.factor(dfResults$MatchID))

par(mfrow=c(length(listMatchIDs),1))

#Plot the AUC

for(strMatchID in listMatchIDs){
  print(strMatchID)
  strMatchID <- as.numeric(as.character(strMatchID))
  pdf(paste0(opt$outdir,"/",opt$db,"-AUC.pdf"))
  dfPlot = dfResults[dfResults$MatchID==strMatchID & dfResults$SyntheticMG=="Illumina_05",]
  print(dim(dfPlot))
  plot(dfPlot$ClustID,dfPlot$MarkerLen,cex=dfPlot$SB_AUC*2,xlab="Clust ID",ylab="Min Marker Length",xlim=c(77,102),main=paste("ShortBRED AUC on the ",opt$db,", Matching at ",strMatchID*100,"% Identity",sep=""),pch=16,col=ifelse(round(dfPlot$SB_AUC,2)==round(max(dfPlot$SB_AUC),2), 2, 1))
  text(dfPlot$ClustID,dfPlot$MarkerLen, round(dfPlot$SB_AUC,2),pos=4,cex=0.8)
  dev.off()
 
  #Plot the Spearman Correlation for ShortBRED
  pdf(paste0(opt$outdir,"/",opt$db,"-corr.pdf"))

  dfPlot = dfResults[dfResults$MatchID==strMatchID & dfResults$SyntheticMG=="Illumina_05",]
  plot(dfPlot$ClustID,dfPlot$MarkerLen,cex=dfPlot$SB_Spearman*2,xlab="Clust ID",ylab="Min Marker Length",xlim=c(77,102),main=paste("ShortBRED Spearman correlation on the ",opt$db,", Matching at ",strMatchID*100,"% Identity",sep=""),pch=16,col=ifelse(round(dfPlot$SB_Spearman,2)==round(max(dfPlot$SB_Spearman),2), 2, 1) )
  text(dfPlot$ClustID, dfPlot$MarkerLen, round(dfPlot$SB_Spearman,2),pos=4,cex=0.8)
  dev.off()
  
  #Plot the Specificity for ShortBRED
  pdf(paste0(opt$outdir,"/",opt$db,"-spec.pdf"))
  
  dfPlot = dfResults[dfResults$MatchID==strMatchID & dfResults$SyntheticMG=="Illumina_05",]
  plot(dfPlot$ClustID,dfPlot$MarkerLen,cex=(dfPlot$SB_Spec)*2,xlab="Clust ID",ylab="Min Marker Length",xlim=c(77,102),main=paste("ShortBRED specificity on the ",opt$db,", Matching at ",strMatchID*100,"% Identity",sep=""),pch=16,col=ifelse(round(dfPlot$SB_Spec,2)==round(max(dfPlot$SB_Spec),2), 2, 1) )
  text(dfPlot$ClustID, dfPlot$MarkerLen, round(dfPlot$SB_Spec,2),pos=4,cex=0.8)
  dev.off()

  #Plot the Sensitivity for ShortBRED
  pdf(paste0(opt$outdir,"/",opt$db,"-sens.pdf"))
  dfPlot = dfResults[dfResults$MatchID==strMatchID & dfResults$SyntheticMG=="Illumina_05",]
  plot(dfPlot$ClustID,dfPlot$MarkerLen,cex=dfPlot$SB_Sens*2,xlab="Clust ID",ylab="Min Marker Length",xlim=c(77,102),main=paste("ShortBRED sensitivity on the ",opt$db,", Matching at ",as.numeric(strMatchID)*100,"% Identity",sep=""),pch=16,col=ifelse(round(dfPlot$SB_Sens,2)==round(max(dfPlot$SB_Sens),2), 2, 1))
  text(dfPlot$ClustID, dfPlot$MarkerLen, round(dfPlot$SB_Sens,2),pos=4,cex=0.8)
  dev.off()
}

