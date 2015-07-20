#! /usr/bin/env Rscript

#Jim Kaminski
rm(list=ls()) 

library(optparse)
#library(stringr)
#library("Biobase")
library(ROCR)
library(ggplot2)

c_ad_ARDB_ROCplot_xlim = c(0,0.2)
c_ad_ARDB_ROCplot_ylim = c(0.80,1.0)

c_ad_VF_ROCplot_xlim = c(0,0.05)
c_ad_VF_ROCplot_ylim = c(0.95,1.0)


c_ad_ARDB_ROCplot_xlim = c(0,0.2)
c_ad_ARDB_ROCplot_ylim = c(0,1)

c_ad_VF_ROCplot_xlim = c(0,0.2)
c_ad_VF_ROCplot_ylim = c(0,1)


c_dSmallConstant <- 10^(-8)
bPctPlot <- FALSE
if(bPctPlot==TRUE){
  c_dSmallConstant <- 0
}

#This variable checks if any value of Predicted or True/Expected
#abundance is lower than our minimal constant. If it's true, then
#we will be missing data on the plot and will want to know.
b_ConstantNotLowest <- FALSE

#####################################################################
#FOR TESTING ONLY, comment this line out when you upload to hutlab3
#setwd("C:\\Users\\Jim\\Desktop\\metagenomes\\")
#####################################################################

#source("http://bioconductor.org/biocLite.R")
#abiocLite("genefilter")

option_list <- list(
  make_option("--A05", default="output/shortbred-sfle/tmp/MG-Multiple-Illumina/ARDB/ARDB05/Results/evaltable.tab", type="character",
              help="ARDB 05%"),
  make_option("--A10", default="output/shortbred-sfle/tmp/MG-Multiple-Illumina/ARDB/ARDB10/Results/evaltable.tab", type="character",
              help="ARDB 10%."),
  make_option("--A25", default="output/shortbred-sfle/tmp/MG-Multiple-Illumina/ARDB/ARDB25/Results/evaltable.tab", type="character",
              help="ARDB 25%"),
  make_option("--V05", default="output/shortbred-sfle/tmp/MG-Multiple-Illumina/VF/VF05/Results/evaltable.tab", type="character",
              help="VFDB 05%"),
  make_option("--V10", default="output/shortbred-sfle/tmp/MG-Multiple-Illumina/VF/VF10/Results/evaltable.tab", type="character",
              help="VFDB 10%"),
  make_option("--V25", default="output/shortbred-sfle/tmp/MG-Multiple-Illumina/VF/VF25/Results/evaltable.tab", type="character",
              help="VFDB 25%"),
  make_option("--cor", default="output/shortbred-sfle/figures/Fig3-AllCorr.png", type="character",
              help=""),
  make_option("--ROC", default="output/shortbred-sfle/figures/Fig2-AllROC.png", type="character",
              help="path for ROC plot"),
  make_option("--AUC", default="output/shortbred-sfle/tmp/MG-Multiple/figures/AUC.tab", type="character",
              help="Output txt file"),
  make_option("--reads", default=5000000, type="numeric",
              help="Enter the number of reads in the metagenome")
  
)


opt <- parse_args(OptionParser(option_list=option_list))
#print(opt$hmp)

A05 <- read.table(opt$A05, header=TRUE, sep="\t", row.names=1,
                  as.is=TRUE)
A10 <- read.table(opt$A10, header=TRUE, sep="\t", row.names=1,
                  as.is=TRUE)
A25 <- read.table(opt$A25, header=TRUE, sep="\t", row.names=1,
                  as.is=TRUE)
V05 <- read.table(opt$V05, header=TRUE, sep="\t", row.names=1,
                  as.is=TRUE)
V10 <- read.table(opt$V10, header=TRUE, sep="\t", row.names=1,
                  as.is=TRUE)
V25 <- read.table(opt$V25, header=TRUE, sep="\t", row.names=1,
                  as.is=TRUE)


###############################################################
#Correlation

strList <- c("A05","A10","A25","V05","V10","V25")
ai_ltype <- c(1,2,3,1,2,3)

astrNames = c("ARDB 5% - 150 genes","ARDB 10% - 500 genes","ARDB 25% - 1000 genes","VFDB 5% - 150 genes","VFDB 10% - 500 genes","VFDB 25% - 1000 genes")
i <- 1




# Make the 
png(file=opt$cor,width=1600,height=1200)
par(mfrow=c(2,3),cex=1.5,mar=(c(5, 4, 6, 2)),las=1,yaxp=c(0.001,0.1,1),ylog=TRUE)




for (dfMG in strList){
  print(dfMG)
  
  dfRes <- get(dfMG)
  
  print(sum(dfRes$WGS.Reads))
  print(sum(dfRes$WGS.Reads)/5000000)
  print(table(dfRes$WGS.Reads>0))
  

  
  print(paste("Plotting...",dfMG,astrNames[i]))
  
  #Abundances
  dfRes$CentAbund <- dfRes$Q.Centroids/sum(dfRes$Q.Centroids)
  dfRes$SBAbund <- dfRes$Q.Markers/sum(dfRes$Q.Markers)
  
  print("ShortBRED")
  MarkerCor <- cor((dfRes$RelCount),dfRes$SBAbund,method="spearman")
  
  print("Centroids")
  CentCor <- cor((dfRes$RelCount),dfRes$CentAbund,method="spearman")
  
  c_dSmallConstant
  
  #Change 0's to a very small constant so so that they display
  
  dfRes$CentAbund[dfRes$CentAbund==0] <- c_dSmallConstant
  dfRes$SBAbund[dfRes$SBAbund==0] <- c_dSmallConstant
  dfRes$RelCount[dfRes$RelCount==0] <- c_dSmallConstant
  
  print(paste("The lowest value really is",min(dfRes$SBAbund,dfRes$CentAbund,dfRes$RelCount))) 
  if (min(dfRes$SBAbund,dfRes$CentAbund,dfRes$RelCount) < c_dSmallConstant){
    b_ConstantNotLowest <- TRUE
  }
  
  dHighPoint <- max(dfRes$RelCount,dfRes$CentAbund,dfRes$SBAbund) 
  print(paste("The high point is ",dHighPoint))
  dHighPoint <- dHighPoint+ .001
  
  

if (bPctPlot==TRUE){
  plot(dfRes$RelCount,dfRes$CentAbund,
       col="blue",main=paste(astrNames[i],opt$name,sep=""),
       pch=3,xlab="Expected Abundance",ylab="Predicted Abundance",ylim=c(c_dSmallConstant,dHighPoint),
       xlim=c(c_dSmallConstant,dHighPoint))
  text(dHighPoint*.6,dHighPoint*.10,paste(" SB :", round(MarkerCor,3)),pos=4,col="red")
  text(dHighPoint*.6,dHighPoint*.05,paste("\n Cen:" ,round(CentCor,3)),pos =4,col="blue" )
}  
else {
  dHighPoint <- 10^(-1)
  plot(dfRes$RelCount,dfRes$CentAbund,
       col="blue",main=paste(astrNames[i],opt$name,sep=""),
       pch=3,xlab="Expected Abundance",ylab="Predicted Abundance",ylim=c(c_dSmallConstant,dHighPoint),
       xlim=c(c_dSmallConstant,dHighPoint),log="xy",yaxt="n",xaxt="n")
  axis(2,at=c(c_dSmallConstant,10^-7,10^-5,10^-3,10^-1),labels=c(0,"1e-07","1e-05","1e-03","1e-01"))
  axis(1,at=c(c_dSmallConstant,10^-7,10^-5,10^-3,10^-1),labels=c(0,"1e-07","1e-05","1e-03","1e-01"))
  text(10^(-3),10^(-7.2),paste(" SB :", round(MarkerCor,3)),pos=4,col="red")
  text(10^(-3),10^(-7.4),paste("\n Cen:" ,round(CentCor,3)),pos =4,col="blue" )
}
  if(i==1){
    legend("topleft",c("ShortBRED","Centroids"),col=c("red","blue"),pch=c(1,3))
  }
  
  points(dfRes$RelCount,(dfRes$SBAbund),col="red",pch=1)
  abline(0,1)
  
  
  i <- i + 1
}
title( "Predicted and Expected Relative Abundance", outer = TRUE,line=2 )
dev.off()

i <- 1

fOut<-file(opt$AUC,open="w")
strHeader <- paste("Metagenome","ShortBRED AUC","Centroids AUC","ShortBRED Sens","Centroid Sens","ShortBRED Spec","Centroid Spec",sep="\t")
writeLines(strHeader,fOut)

ad_SBauc <- c()
ad_Fullauc <- c()


png(filename=opt$ROC,width = 480*2, height = 480)
par(mfrow=c(1,2),las=1)
for (dfMG in strList){
  
  dfRes <- get(dfMG)
  
  #Make labels for ROCR 
  dfRes$PosMarker[dfRes$Q.Markers>0] <- 1
  dfRes$PosMarker[dfRes$Q.Markers==0] <- 0
  
  dfRes$PosCent[dfRes$Q.Centroids>0] <- 1
  dfRes$PosCent[dfRes$Q.Centroids==0] <- 0
  
  dfRes$labels[dfRes$RelCount>0]<- 1
  dfRes$labels[dfRes$RelCount==0]<- 0
  
  predSB <- prediction(dfRes$Q.Markers, dfRes$labels)
  perfSB <- performance(predSB,"tpr","fpr")
  aucSB <- performance(predSB,"auc")
  
  predFull <- prediction(dfRes$Q.Centroids, dfRes$labels)
  perfFull <- performance(predFull,"tpr","fpr")
  aucFull <- performance(predFull,"auc")
  
  #Start new plot when you get the 5% metagenome
  if (dfMG=="A05" | dfMG =="V05"){
    b_add <- FALSE
  }
  else {
    b_add <- TRUE
  }
  
  if(substring(dfMG, 1, 1) == "A"){
    ad_xlim = c_ad_ARDB_ROCplot_xlim
    ad_ylim = c_ad_ARDB_ROCplot_ylim
    str_title = "ARDB"
  }
  else{
    ad_xlim <- c_ad_VF_ROCplot_xlim 
    ad_ylim <- c_ad_VF_ROCplot_ylim 
    str_title = "VF"
  }
  
  
  plot(perfFull,main=str_title,xlim=ad_xlim,ylim=ad_ylim,col="blue",add=b_add,lty=ai_ltype[i],lwd=2)
  plot(perfSB,main=str_title,xlim=ad_xlim,ylim=ad_ylim,col="red",add=TRUE,lty=ai_ltype[i],lwd=2)    
  
  
  # Calculate sensitivity and specificity at 0
  dfRes$TN <- 0
  dfRes$TN[dfRes$labels==0 & dfRes$PosMarker==0] <- 1
  dfRes$FP <- 0
  dfRes$FP[dfRes$labels==0 & dfRes$PosMarker==1] <- 1
  dfRes$TP <- 0
  dfRes$TP[dfRes$labels==1 & dfRes$PosMarker==1] <- 1
  dfRes$FN <- 0
  dfRes$FN[dfRes$labels==1 & dfRes$PosMarker==0] <- 1
  
  dfRes$CentTN <- 0
  dfRes$CentTN[dfRes$labels==0 & dfRes$PosCent==0] <- 1
  dfRes$CentFP <- 0
  dfRes$CentFP[dfRes$labels==0 & dfRes$PosCent==1] <- 1
  dfRes$CentTP <- 0
  dfRes$CentTP[dfRes$labels==1 & dfRes$PosCent==1] <- 1
  dfRes$CentFN <- 0
  dfRes$CentFN[dfRes$labels==1 & dfRes$PosCent==0] <- 1
  
  MarkerSpec <- signif(sum(dfRes$TN) / (sum(dfRes$TN)+sum(dfRes$FP)))
  CentSpec <- signif(sum(dfRes$CentTN) / (sum(dfRes$CentTN)+sum(dfRes$CentFP)))
  
  MarkerSens <- signif(sum(dfRes$TP) / (sum(dfRes$TP)+sum(dfRes$FN)))
  CentSens <- signif(sum(dfRes$CentTP) / (sum(dfRes$CentTP)+sum(dfRes$CentFN)))
  
  strOutput <- paste(astrNames[i],aucSB@y.values,aucFull@y.values,MarkerSens,CentSens,MarkerSpec,CentSpec,sep="\t")
  
  ad_SBauc <- append(ad_SBauc,round(as.numeric(aucSB@y.values),digits=3)*100)
  ad_Fullauc <- append(ad_Fullauc,round(as.numeric(aucFull@y.values),digits=3)*100)
  #aucFull@y.values) 
  
  writeLines(strOutput, fOut)
  
  if(i==3){
    legend("bottomright",c(paste("MG ","SB ","Cent"),paste('5%  ',ad_SBauc[1],ad_Fullauc[1]),
                           paste('10%',ad_SBauc[2],ad_Fullauc[2]),
                           paste('25%',ad_SBauc[3],ad_Fullauc[3])),lty=c(0,1,2,3),lwd=2)
  }
  
  i=i +1
}
close(fOut)
legend("bottomright",c(paste("MG ","SB ","Cent"),paste('5%  ',ad_SBauc[4],ad_Fullauc[4]),
                       paste('10%',ad_SBauc[5],ad_Fullauc[5]),
                       paste('25%',ad_SBauc[6],ad_Fullauc[6])),lty=c(0,1,2,3),lwd=2)
dev.off()


###################################################################################
# Figure 2

par(mfcol=c(2,2))





warnings()
if (b_ConstantNotLowest){
  print("You must set a lower constant. Some points are lower.")
}

