#! /usr/bin/env Rscript

#install.packages("optparse")
rm(list = ls())
library(optparse)
library(stringr)
library(ROCR)

option_list <- list(
  make_option("--name", default="VF", type="character",
              help="Enter the name of the evaluation database (ARDB or VF)."),
  make_option("--rf",  default="C:\\Users\\Jim\\Desktop\\sfletest\\evaltable.txt", type="character",
              help="Enter the name of the results file."),
  make_option("--corr",  default="corr.png", type="character",
              help="Enter the name of the scatterplot png file."),
  make_option("--spec",  default="spec.png", type="character",
              help="Enter the name of the spec png file."),
  make_option("--stats",  default="stats.txt", type="character",
              help="Enter the name of the output stats file."),
  make_option("--tpfp",  default="tprfpr.png", type="character",
              help="Enter the name of the tprfpr png file."),
  make_option("--pct",  default=".99", type="character",
              help="Enter the pct id."),
  make_option("--ml",  default="20", type="character",
              help="Enter the marker length."),
  make_option("--id",  default="99", type="character",
              help="Enter the matching id used."),
  make_option("--metagenome",  default="", type="character",
              help="Enter the name of the metagenome used.")
)

opt <- parse_args(OptionParser(option_list=option_list))

dfResults <- read.table(opt$rf,header=TRUE,sep="\t")


#Make labels for ROCR 
dfResults$PosMarker[dfResults$Q.Markers>0] <- 1
dfResults$PosMarker[dfResults$Q.Markers==0] <- 0

dfResults$PosCent[dfResults$Q.Centroids>0] <- 1
dfResults$PosCent[dfResults$Q.Centroids==0] <- 0

#Plot specificity by cutoff
png(filename=opt$spec,width=600,height=600)
par(mfrow=c(1,1))
dfResults$labels[dfResults$RelCount>0]<- 1
dfResults$labels[dfResults$RelCount==0]<- 0

predSB <- prediction(dfResults$Q.Markers, dfResults$labels)
perfSBspec <- performance(predSB,"spec")
plot(perfSBspec,main="Specificity (presence/absence of family)",xlim=c(0,0.25),col="red")

predFull     <- prediction(dfResults$Q.Centroids, dfResults$labels)
perfFullspec <- performance(predFull,"spec")
plot(perfFullspec,xlim=c(0,0.25),add=TRUE,col="blue")
legend("bottomright",c("ShortBRED","Centroids"),col=c("red","blue"),lty=1)

dev.off()

#Plot TPR by FPR
png(filename=opt$tpfp,height=600,width=600)
par(mfrow=c(1,1))
predSB <- prediction(dfResults$Q.Markers, dfResults$labels)
perfSB <- performance(predSB,"tpr","fpr")
predFull     <- prediction(dfResults$Q.Centroids, dfResults$labels)
perfFull <- performance(predFull,"tpr","fpr")

plot(perfSB,main="ROC Curve",xlim=c(0,0.1), col="red")
plot(perfFull,xlim=c(0,0.1),add=TRUE,col="blue")
legend("bottomright",c("ShortBRED","Centroids"),col=c("red","blue"),lty=1)
dev.off()


#AUC for TPR by FPR
predSB <- prediction(dfResults$Q.Markers, dfResults$labels)
SBauc <- performance(predSB,"auc")
predFull     <- prediction(dfResults$Q.Centroids, dfResults$labels)
Fullauc <- performance(predFull,"auc")

dfResults$TN <- 0
dfResults$TN[dfResults$labels==0 & dfResults$PosMarker==0] <- 1
dfResults$FP <- 0
dfResults$FP[dfResults$labels==0 & dfResults$PosMarker==1] <- 1
dfResults$TP <- 0
dfResults$TP[dfResults$labels==1 & dfResults$PosMarker==1] <- 1
dfResults$FN <- 0
dfResults$FN[dfResults$labels==1 & dfResults$PosMarker==0] <- 1

dfResults$CentTN <- 0
dfResults$CentTN[dfResults$labels==0 & dfResults$PosCent==0] <- 1
dfResults$CentFP <- 0
dfResults$CentFP[dfResults$labels==0 & dfResults$PosCent==1] <- 1
dfResults$CentTP <- 0
dfResults$CentTP[dfResults$labels==1 & dfResults$PosCent==1] <- 1
dfResults$CentFN <- 0
dfResults$CentFN[dfResults$labels==1 & dfResults$PosCent==0] <- 1


MarkerCor <- cor(dfResults$Q.Markers, dfResults$RelCount, method="spearman")
CentCor <- cor(dfResults$Q.Centroids, dfResults$RelCount, method="spearman")
print(MarkerCor)
print(CentCor)

MarkerSpec <- signif(sum(dfResults$TN) / (sum(dfResults$TN)+sum(dfResults$FP)))
CentSpec <- signif(sum(dfResults$CentTN) / (sum(dfResults$CentTN)+sum(dfResults$CentFP)))

MarkerSens <- signif(sum(dfResults$TP) / (sum(dfResults$TP)+sum(dfResults$FN)))
CentSens <- signif(sum(dfResults$CentTP) / (sum(dfResults$CentTP)+sum(dfResults$CentFN)))

fOut<-file(opt$stats)
writeLines(paste(opt$name,opt$metagenome,opt$pct,opt$ml,opt$id,SBauc@y.values,Fullauc@y.values,MarkerCor,CentCor,MarkerSpec, CentSpec,MarkerSens,CentSens,sep="\t"), fOut)
close(fOut)

#sink()
#print(paste(opt$name,"PCTID","MarkerLen",SBauc@y.values,Fullauc@y.values,MarkerCor,CentCor))
#sink()

#print(paste("ShortBRED AUC:",SBauc@y.values))
#print(paste("Centroids AUC:",Fullauc@y.values))


#Values for a cutoff of 0

#print(paste("Specificity at 0 for Markers :"))
#print(signif(sum(dfResults$TN) / (sum(dfResults$TN)+sum(dfResults$FP))),4 )

#print(paste("Specificity at 0 for Centroids :"))
#print(signif(sum(dfResults$CentTN) / (sum(dfResults$CentTN)+sum(dfResults$CentFP))),4 )

#print(paste("Sensitivity at 0 for Markers :"))
#print(signif(sum(dfResults$TP) / (sum(dfResults$TP)+sum(dfResults$FN))),4 )

#print(paste("Sensitivity at 0 for Centroids :"))
#print(signif(sum(dfResults$CentTP) / (sum(dfResults$CentTP)+sum(dfResults$CentFN))),4 )

#table(dfResults$PosMarker)
#table(dfResults$PosMarker,dfResults$labels)

#print(paste("Correlation of Marker Abundance with True Count :"))

#print(MarkerCor)

#print(paste("Correlation of Centroid Abundance with True Count :"))

#print(CentCor)




attach(dfResults)
dHighPoint <- max(RelCount,Q.Centroids/sum(Q.Centroids),Q.Markers/sum(Q.Markers)) 
png(filename=opt$corr,width=600,height=600)
plot(RelCount,Q.Markers/sum(Q.Markers),
     col="red",main=paste("Predicted and Actual Relative Abundance: ",opt$name,sep=""),
     pch=20,xlim=c(0,dHighPoint),ylim=c(0,dHighPoint),xlab="Actual",ylab="Predicted")
legend("topleft",c("ShortBRED","Centroids"),col=c("red","blue"),pch=15)
text(0.04,0.03,paste("Mar:", round(MarkerCor, 3), "\n Cen:" ,round(CentCor,3) ),pos =4 )
points(RelCount,Q.Centroids/sum(Q.Centroids),col="blue",pch=20)
text(RelCount[High==1],SBRelCount[High==1],Family[High==1], cex=0.8, pos=4, col="black")
abline(0,1)

dev.off()
