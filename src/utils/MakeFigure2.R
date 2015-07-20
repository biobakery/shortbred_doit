#! /usr/bin/env Rscript

#Jim Kaminski
rm(list=ls()) 

library(optparse)
#install.packages("gridExtra")
#install.packages("reshape",lib="/n/home11/jkaminski/Rpackages")
library(gridExtra,lib.loc="/n/home11/jkaminski/Rpackages")
#library(stringr)
#library("Biobase")
library(ROCR)
library(ggplot2)
library(reshape,lib.loc="/n/home11/jkaminski/Rpackages")



theme_set(theme_gray(base_size = 8))

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
  make_option("--fig", default="output/shortbred-sfle/figures/Figure2-Illumina.pdf", type="character",
              help=""),
  make_option("--AUC", default="output/shortbred-sfle/tmp/MG-Multiple/figures/AUC.tab", type="character",
              help="Output txt file"),
  make_option("--reads", default=5000000, type="numeric",
              help="Enter the number of reads in the metagenome"),
  make_option("--y_lowARDB", default=0.80, type="numeric",
              help="Enter the number of reads in the metagenome"),
  make_option("--y_lowVFDB", default=0.95, type="numeric",
              help="Enter the number of reads in the metagenome"),
  make_option("--title", default="ShortBRED", type="character",
              help="Enter the title of the figure."),
  make_option("--fullROC", default="N", type="character",
              help="Enter Y if you want full axes on ROC plot."),
  make_option("--eps", default="fig.eps", type="character",
              help="Enter the name for the EPS file.")
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


strMainTitle <- opt$title
				  
fOut<-file(opt$AUC,open="w")
strHeader <- paste("Metagenome","ShortBRED AUC","Centroids AUC","ShortBRED Sens","Centroid Sens","ShortBRED Spec","Centroid Spec",
                   "ShortBRED Corr","Centroid Corr",sep="\t")
writeLines(strHeader,fOut)

ad_SBauc <- c()
ad_Fullauc <- c()

if ( opt$fullROC == "Y" ){
  print("Making large ROC")
	c_ad_ARDB_ROCplot_xlim = c(0,0.2)
	c_ad_ARDB_ROCplot_ylim = c(0.5,1.0)

	c_ad_VF_ROCplot_xlim = c(0,0.2)
	c_ad_VF_ROCplot_ylim = c(0.5,1.0)
} else{
	c_ad_ARDB_ROCplot_xlim = c(0,0.2)
	c_ad_ARDB_ROCplot_ylim = c(.8,1.0)

	c_ad_VF_ROCplot_xlim = c(0,0.05)
	c_ad_VF_ROCplot_ylim = c(0.95,1.0)
}





#################################################################
# Make the ROC Plots

ai_ltype <-c(1,2,3)
astrNames<-c("5% - 150 genes","10% - 500 genes", "25% - 1000 genes")

ROCPlot <- function(strList,strTitle,ad_xlim,ad_ylim,fOut,astrNames){
  
  ggplot_ROC <- ggplot()
  i<- 0

  dfAllROC <- data.frame( Meth=character(),MG=character(),
                          MethAndMG= character(), x = numeric(), y = numeric())
  print(dfAllROC)
  for (dfMG in strList){
    i <- i +1
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
    
    dfRes$CentAbund <- dfRes$Q.Centroids/sum(dfRes$Q.Centroids)
    dfRes$SBAbund <- dfRes$Q.Markers/sum(dfRes$Q.Markers)
    
    MarkerCor <- cor((dfRes$RelCount),dfRes$SBAbund,method="spearman")
    CentCor <- cor((dfRes$RelCount),dfRes$CentAbund,method="spearman")
    
    strOutput <- paste(strList[i],aucSB@y.values,aucFull@y.values,MarkerSens,CentSens,MarkerSpec,CentSpec,MarkerCor,CentCor,sep="\t")
    
    ad_SBauc <- append(ad_SBauc,round(as.numeric(aucSB@y.values),digits=3)*100)
    ad_Fullauc <- append(ad_Fullauc,round(as.numeric(aucFull@y.values),digits=3)*100)
    #aucFull@y.values) 
    
    writeLines(strOutput, fOut)
    
    #Convert peformance objects into dataframes. This makes it easier to use ggplot2
    dfPerfFull <- data.frame(Meth="Centroids",
                             MG=paste(substr(dfMG,2,3),"%",sep=""),
                            MethAndMG=paste(dfMG,"_Cent",sep=""),
                             x=perfFull@x.values[[1]],y=perfFull@y.values[[1]])
    dfPerfSB <- data.frame(Meth="ShortBRED",MG=paste(substr(dfMG,2,3),"%",sep=""),MethAndMG=paste(dfMG,"_SB",sep=""),x=perfSB@x.values[[1]],y=perfSB@y.values[[1]])  
    
    dfAllROC<-  rbind(dfAllROC,dfPerfSB)
    dfAllROC<-  rbind(dfAllROC,dfPerfFull)
    print("This ran")
  }
    print(dfAllROC)
    ggplot_ROC <- ggplot(data=dfAllROC,aes(x=x,y=y,linetype=MG,colour=Meth),size=1,lwd=1) +
    scale_colour_manual(values=c("red","blue"),guide = guide_legend(title = "Method")) +
    scale_linetype_manual(values = c(1,2,3),
          guide = guide_legend(title =paste("Synthetic","metagenome","proportion",sep='\n'))) +
    geom_line() +
    coord_cartesian(xlim=ad_xlim, ylim=ad_ylim) +  
    xlab("False Positive Rate") +
    ylab("True Positive Rate") + labs(title = strTitle) +
    theme(legend.key.size= unit(0.25, "cm"))
  
 
  
    return(ggplot_ROC)
}
  


#################################################################
# Make the Correlation Plots

CorrPlot <- function(dfRes,strTitle){
  c_dSmallConstant <- 10^(-8)
  
      
  print(paste("Plotting..."))
    
  #Abundances
  dfRes$CentAbund <- dfRes$Q.Centroids/sum(dfRes$Q.Centroids)
  dfRes$SBAbund <- dfRes$Q.Markers/sum(dfRes$Q.Markers)
    
  print("ShortBRED")
  MarkerCor <- cor((dfRes$RelCount),dfRes$SBAbund,method="spearman")
  #print(MarkerCor)  
  print("Centroids")
  CentCor <- cor((dfRes$RelCount),dfRes$CentAbund,method="spearman")
  
  #Take the "true" abundance,SB abundance and ShortBRED abundance, 
  dfPlotResults <-dfRes[,c("RelCount","CentAbund","SBAbund")]
  colnames(dfPlotResults) <- c("RelCount","Centroids","ShortBRED")
  dfPlotResults <- melt(dfPlotResults,id=c("RelCount"), measured=c("ShortBRED","Centroids"))
  colnames(dfPlotResults)[colnames(dfPlotResults)=="variable"] <- "Method"
  #print(dfPlotResults)
  
  strSpearman = paste("ShortBRED: ",round(MarkerCor,digits=3),"\n","Centroids: ",round(CentCor,digits=3),sep="")
  
  dHighPoint <- max((dfRes$CentAbund),(dfRes$SBAbund))+.0001
  ggplot_Corr<- ggplot(data=dfPlotResults,aes(x=RelCount,y=value,colour=Method,shape=Method)) + geom_point() +
    scale_colour_manual(values=c("blue","red")) + scale_shape_manual(values=c(3,1)) +
    scale_y_sqrt() + scale_x_sqrt() +
    coord_cartesian(xlim=c(0.00,dHighPoint), ylim=c(0.00,dHighPoint)) +
    xlab("Predicted Abundance") +
    ylab("Expected Abundance") + labs(title = strTitle) +
    geom_abline() + theme(legend.key.size= unit(0.25, "cm")) +
    annotate("text", label = strSpearman, x = dHighPoint/2, y = dHighPoint/30, size = 2, colour = "black")
 
  return(ggplot_Corr)
}


ggplot_ARDB_ROC <- ROCPlot(strList=c("A05","A10","A25"),strTitle="A. Antibiotic Resistance Genes Database \n ROC",ad_xlim=c_ad_ARDB_ROCplot_xlim,
                       ad_ylim=c_ad_ARDB_ROCplot_ylim,fOut=fOut)
ggplot_VFDB_ROC <- ROCPlot(strList=c("V05","V10","V25"),strTitle="B. Virulence Factors Database \n ROC",ad_xlim=c_ad_VF_ROCplot_xlim,
                       ad_ylim=c_ad_VF_ROCplot_ylim,fOut=fOut)

ggplot_ARDB_CORR <- CorrPlot(A10,strTitle="C. Antibiotic Resistance Genes Database \n Correlation - 10% of Metagenome, 500 genes ")
ggsave(file="test.pdf",ggplot_ARDB_CORR)
ggplot_VFDB_CORR <- CorrPlot(V10,strTitle="D. Virulence Factors Database \n Correlation - 10% of Metagenome, 500 genes ")

grid_Fig2 <- arrangeGrob(ggplot_ARDB_ROC,ggplot_VFDB_ROC, ggplot_ARDB_CORR,
                         ggplot_VFDB_CORR, nrow=2)

#grid.draw(grid_Fig2)
#print(class(grid_Fig2))
ggsave(file=opt$fig,grid_Fig2,width=9)
ggsave(file=opt$eps,grid_Fig2,width=7.0,height=6,units="in")

