#! /usr/bin/env Rscript


#Jim Kaminski
rm(list=ls()) 

print("START")

#source("http://bioconductor.org/biocLite.R") 
#biocLite("Biobase")
library(Biobase)
library(methods)
require(Biobase)

#install.packages("RColorBrewer")
#library(RColorBrewer)
library(optparse)
library(stringr)
library("gplots")

getAnywhere("new")
new("ExpressionSet")
sessionInfo(package = "Biobase")
####################################################################
#Constants

file.sep <- .Platform$file.sep

####################################################################
#FOR TESTING ONLY, comment this line out when you upload to hutlab3
#setwd("C:\\Users\\Jim\\Desktop\\ProcessSB\\")
#####################################################################

#source("http://bioconductor.org/biocLite.R")

option_list <- list(
  make_option("--hmp", default="hmp_merged\\merged.txt", type="character",
              help="Enter the HMP SB results file."),
  make_option("--hmet", default="Export\\hmp_metadata.dat", type="character",
              help="Enter the HMP metadata."),
  make_option("--hreads", default="summary_nreads.txt", type="character",
              help="Enter the HMP reads file."),
  make_option("--yat", default="yat_merged\\merged.txt", type="character",
              help="Enter the Yatsunenko SB results file."),
  make_option("--ymet", default="YatMap.txt", type="character",
              help="Enter the Yatsunenko metadata file."),
  make_option("--ardb", default="ARDB.txt", type="character",
              help="Enter the ARDB."),
  make_option("--out", default="C:\\Users\\Jim\\Desktop\\ProcessSB\\", type="character",
              help="output directory.")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt$hmp)

#################################################################################
#MakeARDB() - Function to make dfARDBMeta

MakeARDB <- function(){
  dfGenes <- read.table(paste(opt$ardb,file.sep,"tabs",file.sep,"abrg.tab",sep=""))
  dfGenes <- dfGenes[,c(1:2)]
  colnames(dfGenes) <- c("GeneName","Type")
  
  
  dfOrigin <- read.table(paste(opt$ardb,file.sep,"tabs",file.sep,"originType.tab",sep=""))
  dfOrigin <- dfOrigin[,c(1,3)]
  colnames(dfOrigin) <-c("Type","Class")
  
  dfClassInfo <- read.table(paste(opt$ardb,file.sep,"tabs",file.sep,"classInfo.tab",sep=""),sep="\t",allowEscapes=TRUE)
  dfClassInfo[,2]<- str_replace(dfClassInfo[,2],"\n","")
  #write.table(dfClassInfo, file="Export\\CLASSES.txt",sep="\t")
  colnames(dfClassInfo) <-c ("Class","Description")
  
  
  dfARDBmeta <- merge(dfGenes,dfOrigin,by="Type",all.x=TRUE,sort=TRUE)                    
  dfARDBmeta$Class <- tolower(dfARDBmeta$Class)
  dfARDBmeta <- merge(dfARDBmeta,dfClassInfo,by="Class",all.x=TRUE,sort=TRUE)
  
  rownames(dfARDBmeta)<-dfARDBmeta$GeneName
  
  print("ARDB Loaded....")
  return(dfARDBmeta)
  
}

##################################################################################
#Prepare ARDB metadata
dfARDBmeta <- MakeARDB()

#################################################################
#Load in the HMP data
matHMP <- as.matrix(read.table(opt$hmp, header=TRUE, sep="\t", row.names=1,
                               as.is=TRUE))


matHMP[(is.na(matHMP))] <- 0
print("HMP Results Loaded....")

#################################################################
#Keep anything ending in "0", that indicates the ShortBRED counts

print(str_match(colnames(matHMP),".*0$")[,1])
matHMP <- matHMP[,colnames(matHMP) %in% str_match(colnames(matHMP),".*0$")[,1]]

print(colnames(matHMP))

#############################################################################
#Make HMP ExpressionSet

#Process HMP Metadata
dfHMPMeta <- as.data.frame((read.table(opt$hmet, row.names=1,header=TRUE, sep="\t")))
iCols = ncol(dfHMPMeta)
print("Loaded HMP Metadata...")

#Cut down the metadata to the samples we profiled.
#Note that this code is currently set to just get NReads for the HMP
dfHMPMeta <- merge(dfHMPMeta,t(matHMP),by="row.names",all.y=TRUE)
rownames(dfHMPMeta) <- dfHMPMeta[,1]
dfHMPMeta <- dfHMPMeta[,c(1:iCols)]
names(dfHMPMeta)[names(dfHMPMeta)=="Row.names"] <- "SN"

dfHMPReads <- as.data.frame((read.table(opt$hreads, header=FALSE, sep="\t")))
colnames(dfHMPReads) <- c("SN","NReads")
dfHMPMeta <- merge(dfHMPMeta,dfHMPReads,by="SN",all.x=TRUE)

#Confirm that the set of samples is the same for both.
#The ExpressionSet function demands it.
rownames(dfHMPMeta)=dfHMPMeta$SN
print(rownames(dfHMPMeta))
print(colnames(matHMP))
print(all(rownames(dfHMPMeta)==colnames(matHMP)))

#Cut ARDB file down to genes found in HMP. ExpressionSet demands it.
rownames(dfARDBmeta)<-dfARDBmeta$GeneName
dfARDBmeta <- dfARDBmeta[(dfARDBmeta$GeneName %in% rownames(matHMP)) ,]
dfARDBmeta <- dfARDBmeta[order(dfARDBmeta$GeneName) , ]

#Make the expression set.
#expHMP <- ExpressionSet(assayData=as.matrix(matHMP),
#                        phenoData=new("AnnotatedDataFrame",data=dfHMPMeta),
#                        featureData =  new("AnnotatedDataFrame",data=dfARDBmeta))

#print ("Made HMP Expression Set....")
################################################################################
#Export HMP Expression Set

#ARDB + Counts
#dfCountsClass<- merge(dfARDBmeta,dfHMPsummed,by="row.names")
#write.table(dfCountsClass, file="Export\\HMP_SB.txt",sep="\t",row.names=FALSE)

#MetaData
#dfHMP_TMeta <- as.data.frame(t(dfHMPMeta))
#write.table(dfHMP_TMeta, file="Export\\HMP_Meta.txt",sep="\t",col.names=FALSE)




#############################################################################
#Make YAT ExpressionSet
print("Working on Yatsunenko data....")
#Load in the YAT data
matYAT <- as.matrix(read.table(opt$yat, header=TRUE, sep="\t", row.names=1,
                               as.is=TRUE))

matYAT[(is.na(matYAT))] <- 0

#Keep columns ending in 0, they contain the ShortBRED counts
print(str_match(colnames(matYAT),".*0$")[,1])
matYAT <- matYAT[,colnames(matYAT) %in% str_match(colnames(matYAT),".*0$")[,1]]

#Rename the columns to just have the sample names
colnames(matYAT) <- str_match(colnames(matYAT), "([0-9]*).fna.*")[,2]


#This observation does not have metadata, so I am dropping it.
#matYAT <- matYAT[,!names(matYAT) %in% "20716"]


#Load in the metadata
dfYATmeta <- as.data.frame(read.table(opt$ymet, row.names=1, header=FALSE, sep="\t"))
iCols = ncol(dfYATmeta)

#Assign column names
colnames(dfYATmeta) <- c("id1","id2","id3","Country","City","Age","Reads","QualReads")

#Cut down metadata to samples that I profiled. Merge on YAT results, then cut down to metadata.
dfYATmeta <- merge(dfYATmeta,t(matYAT),by="row.names",all.y=TRUE)
dfYATmeta <- dfYATmeta[,c(1:iCols)]
rownames(dfYATmeta)<- dfYATmeta[,1]

print("Yatsunenko data loaded..")
##################################################################################
#Prepare ARDB metadata
dfARDBmeta <- MakeARDB()

#Cut down ARDB to genes found in Yatsunenko data.
dfARDBmeta <- dfARDBmeta[(dfARDBmeta$GeneName %in% rownames(matYAT)) ,]
dfARDBmeta <- dfARDBmeta[order(dfARDBmeta$GeneName) , ]
all(rownames(dfARDBmeta) == dfARDBmeta$GeneName)

expYAT <- ExpressionSet(assayData=matYAT,
                        phenoData=new("AnnotatedDataFrame",data=dfYATmeta),
                        featureData =  new("AnnotatedDataFrame",data=dfARDBmeta))

################################################################################
#Export YAT Expression Set

#ARDB + Counts
dfCountsClass<- merge(dfARDBmeta,matYAT,by="row.names")
write.table(dfCountsClass, file="Export\\YAT_SB.txt",sep="\t",row.names=FALSE)

#MetaData
write.table(t(dfYATmeta), file="Export\\YAT_Meta.txt",sep="\t",col.names=FALSE)



###############################################################################
#Normalize ShortBRED values for heatmaps

#HMP
exprs(expHMP)[is.na(exprs(expHMP))]<- 0
exprs(expHMP) <- t(apply(exprs(expHMP),1,"/",c(as.numeric(as.character(phenoData(expHMP)$NReads)))) )
exprs(expHMP)[is.na(exprs(expHMP))]<- 0

#Yatsunenko
exprs(expYAT)[is.na(exprs(expYAT))]<- 0
exprs(expYAT)<- t(apply(exprs(expYAT),1,"/",c(as.numeric(phenoData(expYAT)$Reads))) )
exprs(expYAT)[is.na(exprs(expYAT))]<- 0

###########################################################
#Sum by antibiotic group
matHMP_ABR <- esApply(expHMP, 2, function(x){
  tapply(x, featureData(expHMP)$Class, FUN=sum)
})

matYAT_ABR <- esApply(expYAT, 2, function(x){
  tapply(x, featureData(expYAT)$Class, FUN=sum)
})

#Make these into expression files as well
expYAT_ABR <- ExpressionSet(assayData=matYAT_ABR,
                        phenoData=new("AnnotatedDataFrame",data=dfYATmeta),
                       )

expHMP_ABR <- ExpressionSet(assayData=matHMP_ABR,
                            phenoData=new("AnnotatedDataFrame",data=dfHMPMeta),
)

#I remove observation 20716 because it does not have any metadata.
expYAT <- expYAT[,!sampleNames(expYAT) %in% "20716"]
expYAT_ABR <- expYAT_ABR[,!sampleNames(expYAT_ABR) %in% "20716"]

heatmap(matHMP_ABR,main="HMP",labRow=rownames(matHMP_ABR),scale="col",cexRow=1,cexCol=1,col=heat.colors(256))
heatmap(matYAT_ABR,main="YAT",labRow=rownames(matYAT_ABR),scale="col",cexRow=1,cexCol=1,col=heat.colors(256))

#############################################################################
#Heatmaps 
lstrConds <- c("col","row","none")
hmcols<-colorRampPalette(c("white","red"))(256)

png(paste(opt$out,"HMP.png",sep=""),1000,1000)
heatmap.2(exprs(expHMP_ABR),main="HMP",sub=paste("by","col"),labRow=featureData(expHMP)$Class,
        scale="col",col=hmcols,labCol="",cexRow=1,cexCol=1,trace="none",density.info="none",key = TRUE)
dev.off()


for(cond in lstrConds){
  print(cond)
  width<-1000
  height<-1000
  
  png(paste(opt$out,"HMP-",cond,".png",sep=""),width,height)
  heatmap.2(exprs(expHMP),main="HMP",sub=paste("by",cond),labRow=featureData(expHMP)$Class,
          scale=cond,cexRow=1,cexCol=1,col=hmcols,labCol="",trace="none",density.info="none",key=TRUE)
  dev.off()
  
  png(paste(opt$out,"YAT-",cond,".png",sep=""),width,height)
  heatmap.2(exprs(expYAT),main="YAT",sub=paste("by",cond),labRow=featureData(expYAT)$Class,
          scale=cond,cexRow=1,cexCol=1,col=hmcols,labCol="",trace="none",density.info="none",key=TRUE)
  dev.off()
  
  png(paste(opt$out,"HMPbyClass-",cond,".png",sep=""),width,height)
  heatmap.2(exprs(expHMP_ABR),main="HMP by ABR Class",labRow=rownames(matHMP_ABR),
          scale=cond,cexRow=1,cexCol=1,col=hmcols,labCol="",trace="none",density.info="none",key=TRUE)
  dev.off()
  
  png(paste(opt$out,"YATbyClass-",cond,".png",sep=""),width,height)
  heatmap.2(exprs(expYAT_ABR),main="YAT by ABR Class",labRow=rownames(matYAT_ABR),
          scale=cond,cexRow=1,cexCol=1,col=hmcols,labCol="",trace="none",density.info="none",key=TRUE)
  dev.off()
}

###########################################################
#Remaining Tasks



#featureData(expHMP)$Class
by(exprs(expHMP),featureData(expHMP)$Class,sum)

heatmap(matHMP_ABR,main="HMP by ABR Class",labRow=rownames(matHMP_ABR),
        scale="col",cexRow=1,cexCol=1,col=hmcols,labCol="")


#write.table(levels(factor(featureData(expHMP)$Class)),"ARDBclassesOut.txt", col.names=NA,quote=FALSE,sep="\t")




HMPclasses <- levels(factor(featureData(expHMP)$Class)) 
YATclasses <- levels(factor(featureData(expYAT)$Class))

#MERGE ARDB CLASSES ONTO THE UNIION OF HMP CLASSES AND YATCLASSES
DataClasses <- as.data.frame(union(HMPclasses,YATclasses))

ARDBClasses <- (read.table("TableForARDBDescriptions.txt", header=TRUE, sep="\t",
                           as.is=TRUE))

colnames(DataClasses) <- "Class"
ARDBMergedClasses <-merge(DataClasses,ARDBClasses,by="Class",all.x=TRUE)

write.table(ARDBMergedClasses,"Table2.txt", col.names=NA,quote=FALSE,sep="\t")

YATclasses[!YATclasses %in% HMPclasses]

mean(dfHMPMeta$NReads)


#Check Age
HMPAge <-as.data.frame(sampleNames(expHMP))
rownames(HMPAge) <- HMPAge[,1]

HMPMeta2 <- as.data.frame(t(read.table("HMPreduced.txt", header=TRUE, sep="\t",row.names=1)))

HMPAge <- merge(HMPAge,HMPMeta2,by="row.names")
mean(as.numeric(as.character(HMPAge$AGEENR)))


######################################################
#Wilcoxon Tests

#Test if median Tet_rpp > 0
wilcox.test(matHMP_ABR["tet_rpp",],alternative="greater",conf.int=TRUE)
wilcox.test(matYAT_ABR["tet_rpp",],alternative="greater",conf.int=TRUE)

#Test if median Tet_rpp > median (median positive expression value)
PosHMPMedians <-apply(matHMP_ABR, 2, function(x){median(x[x > 0])})
PosYATMedians <-apply(matYAT_ABR, 2, function(x){median(x[x > 0])})


wilcox.test(matHMP_ABR["tet_rpp",],PosHMPMedians,paired=TRUE,alternative="g",conf.int=TRUE)
wilcox.test(matYAT_ABR["tet_rpp",],PosYATMedians,paired=TRUE,alternative="g",conf.int=TRUE)

#Test if median bla_a > 0
wilcox.test(matHMP_ABR["bla_a",],alternative="greater",conf.int=TRUE)
wilcox.test(matYAT_ABR["bla_a",],alternative="greater",conf.int=TRUE)

wilcox.test(matHMP_ABR["bla_a",],PosHMPMedians,paired=TRUE,alternative="g",conf.int=TRUE)
wilcox.test(matYAT_ABR["bla_a",],PosYATMedians,paired=TRUE,alternative="g",conf.int=TRUE)

#Test if median bla_a > 0
wilcox.test(exprs(matYAT_ABR)["bla_a",],alternative="greater",conf.int=TRUE)

#featureData(expYAT_ABR)$Class=="tet_rpp"

wilcox.test(matYAT_ABR["bla_a",],PosYATMedians,paired=TRUE,alternative="g",conf.int=TRUE)
wilcox.test( exprs(expYAT_ABR[,phenoData(expYAT_ABR)$Country!="USA"]),exprs(expYAT_ABR[,phenoData(expYAT_ABR)$Country=="USA"]),alternative="g",conf.int=TRUE)



######################################3
#Testing

#Yatsunenko
exprs(expYAT)["AAB03644","20717"]*120104

exprs(expYAT)["AAB03644","20751"]*180645
exprs(expYAT)["AAB51122","20751"]*180645 

exprs(expYAT)["P30899","20826"]*381196
exprs(expYAT)["P36890","20826"]*381196
 

#HMP
exprs(expHMP)["AAB03644","SRS011061"]*53514460

exprs(expHMP)["AAB03644","SRS011586"]*115363914
exprs(expHMP)["AAB51122","SRS011586"]*115363914





#testDF <- as(expHMP,"data.frame")
#write.exprs(ExpressionSet)
#write.exprs(ExpressionSet)
#write.exprs(ExpressionSet)