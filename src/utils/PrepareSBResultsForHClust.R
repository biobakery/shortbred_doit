#! /usr/bin/env Rscript
rm(list=ls()) 
library(optparse)



option_list <- list(
  make_option("--all_results",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-05-28\\AllResults+ABRClass.tab", type="character",
              help="Enter the file with detailed ABR results."),
  make_option("--T2D",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\data\\T2D\\T2Dmeta.txt", type="character",
              help="Enter the file with the T2D metadata."),
  make_option("--HMP",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\data\\HMP\\HMPmeta.txt", type="character",
              help="Enter the file with the Yatsunenko metadata."),
  make_option("--YAT",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\data\\YAT\\YatMap.txt", type="character",
              help="Enter the file with the Yatsunenko metadata."),
  make_option("--YATS2",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\data\\YAT\\YatTableS2.txt", type="character",
              help="Enter the file with the Yatsunenko metadata."),
#############################
# Output
#############################
  make_option("--out_results",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-05-28\\DetailedClean.tab", type="character",
              help="Enter the file with the Yatsunenko metadata."),
  make_option("--out_classresults",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-05-28\\ABRClassesClean.tab", type="character",
              help="Enter the file with the Yatsunenko metadata."),
  make_option("--out_meta",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-05-28\\MetaClean.tab", type="character",
              help="Enter the file with the Yatsunenko metadata."),
  make_option("--out_meta_gender",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-05-28\\MetaCleanGender.tab", type="character",
            help="Enter the file with the Yatsunenko metadata."),
  make_option("--out_meta_ds",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-05-28\\MetaCleanDS.tab", type="character",
            help="Enter the file with the Yatsunenko metadata."),
make_option("--out_meta_country",default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-05-28\\MetaCleanCountry.tab", type="character",
            help="Enter the file with the Yatsunenko metadata.")

)
opt <- parse_args(OptionParser(option_list=option_list))

astrClasses <- c("character","character",rep("numeric",552))
dfDetailed <-read.table(opt$all_results,row.names=1,quote="",comment.char="",header=TRUE,sep="\t",colClasses=astrClasses)
plot(colSums(dfDetailed[,-1]))

colnames(dfDetailed) <- sub("AllResults..ARDBfilter_WashU_","",colnames(dfDetailed) )
colnames(dfDetailed) <- sub("^.*results..","",colnames(dfDetailed) )
colnames(dfDetailed) <- sub(".fna..Count","",colnames(dfDetailed),fixed=T )
colnames(dfDetailed) <- sub("RPKM..Count","",colnames(dfDetailed),fixed=T )
colnames(dfDetailed) <- sub("SRS017103.1","SRS017103",colnames(dfDetailed),fixed=T )
colnames(dfDetailed)[1] <- "ABR_Class"
dfClasses <- aggregate(dfDetailed[,-1], by=list(dfDetailed[,1]),FUN=sum, na.rm=TRUE)
colnames(dfClasses)[1]<- "ABR_Class"
colnames(dfClasses) <- sub("SRS017103.1","SRS017103",colnames(dfClasses),fixed=T )
###################################################################################################
# Prepare T2D metadata
###################################################################################################
dfT2DMeta <- read.table(opt$T2D,quote="",comment.char="",header=TRUE,sep="\t")
dfT2DMeta <- dfT2DMeta[,c(1,2,7)]
dfT2DMeta$Country <- "China"
dfT2DMeta$Dataset <- "T2D-N"
dfT2DMeta[dfT2DMeta$Diabetic..Y.or.N.=="Y","Dataset"] <- "T2D-Y"
dfT2DMeta <- dfT2DMeta[,-3]

###################################################################################################
# Prepare YAT metadata
###################################################################################################
dfYATMetaGender <- read.table(opt$YATS2,quote="",comment.char="",header=TRUE,sep="\t")
dfYATMetaGender <-dfYATMetaGender[,c(5,8)]
colnames(dfYATMetaGender) <- c("Meta.ID","Gender")

dfYATMeta <- read.table(opt$YAT,quote="",comment.char="",header=TRUE,sep="\t")
dfYATMeta <- dfYATMeta[,c(1,4,5)]
colnames(dfYATMeta) <- c("Sample.ID","Meta.ID","Country")

dfYATMeta <- merge(dfYATMeta,dfYATMetaGender,by="Meta.ID")

colnames(dfClasses)
dfYATMeta$Dataset <- "Yatsunenko"

#Order Should be (Sample.ID  Gender	Country	Dataset)
dfYATMeta <- dfYATMeta[c("Sample.ID","Gender","Country","Dataset")]
dfYATMeta$Sample.ID <- as.character(dfYATMeta$Sample.ID)
levels(dfYATMeta$Country)[levels(dfYATMeta$Country)=="USA"] <- "USA-Global"

###################################################################################################
# Prepare HMP metadata
###################################################################################################
dfHMPMeta <- read.table(opt$HMP,quote="",comment.char="",header=TRUE,sep="\t")
dfHMPMeta <- dfHMPMeta[,c("SRS","Gender")]
dfHMPMeta$Country <- "USA-HMP"
dfHMPMeta$Dataset <- "HMP"
colnames(dfHMPMeta)<- c("Sample.ID","Gender","Country","Dataset")
###################################################################################################
# Combine the metadata for the three datasets
###################################################################################################

dfMeta <- rbind(dfT2DMeta,dfYATMeta)
dfMeta <- rbind(dfMeta,dfHMPMeta)

###################################################################################################
# Keep only the ID's that are in dfClasses 
###################################################################################################
dfMeta$Sample.ID  <- sub("-",".",dfMeta$Sample.ID ,fixed=T )

table(dfMeta[dfMeta$Sample.ID %in% colnames(dfClasses),"Dataset"])

astrMissingData <- colnames(dfClasses)[!colnames(dfClasses) %in% dfMeta$Sample.ID]
astrMissingData <- astrMissingData[!astrMissingData %in% "ABR_Class"]

for (strID in astrMissingData){
  if (any(grep("CON",strID),grep("T2D",strID))){
    dfRow <- data.frame(Sample.ID=strID,Gender="Unknown",Country="China",Dataset="T2D")
  }
  else{
    dfRow <- data.frame(Sample.ID=strID,Gender="Unknown",Country="Unknown",Dataset="Yatsunenko")
  }
  dfMeta <- rbind(dfMeta,dfRow)
}

dfMeta <- dfMeta[dfMeta$Sample.ID %in% colnames(dfDetailed),]
dfMeta$Sample.ID[duplicated(dfMeta$Sample.ID)]
colnames(dfClasses)[!colnames(dfClasses) %in% dfMeta$Sample.ID]
colnames(dfClasses)[duplicated(colnames(dfClasses))]

levels(dfMeta$Gender)[levels(dfMeta$Gender)=="female"] <- "Female"
levels(dfMeta$Gender)[levels(dfMeta$Gender)=="male"] <- "Male"
levels(dfMeta$Gender)[levels(dfMeta$Gender)=="female "] <- "Female"
levels(dfMeta$Gender)[levels(dfMeta$Gender)=="male "] <- "Male"
levels(dfMeta$Gender)[levels(dfMeta$Gender)=="NA"] <- "Unknown"
table(dfMeta$Dataset)
table(dfMeta$Country)
# Write out tables with only gender metadata, or only Dataset metadata
dfGender <- dfMeta[,c("Sample.ID","Gender")]
dfDS <- dfMeta[,c("Sample.ID","Dataset")]
dfCountry <- dfMeta[,c("Sample.ID","Country")]
write.table(dfGender,opt$out_meta_gender,sep="\t",row.names=F,col.names=F,quote=F)
write.table(dfDS,opt$out_meta_ds,sep="\t",row.names=F,col.names=F,quote=F)
write.table(dfCountry,opt$out_meta_country,sep="\t",row.names=F,col.names=F,quote=F)

###################################################################################################
# Transpose metadata
##############################################################################################
rownames(dfMeta)<-dfMeta$Sample.ID

AddMetaData <- function(dfData,dfMeta){
  matData <- t(dfData)
  dfMerged <- merge(matData,dfMeta,by="row.names",all.x=T)
  dfMerged <- as.data.frame(t(dfMerged))
  #colnames(dfMerged) <- dfMerged$Row.names
  return(dfMerged)
  
}
dfDetClean <- AddMetaData(dfDetailed,dfMeta)
dfClassClean <- AddMetaData(dfClasses,dfMeta)


##################################################################################################
# Additional Cleanup for Classes File
##################################################################################################

iIndex <- which(dfClassClean =="ABR_Class", arr.ind = T)[2]
aiColOrder <- 1:dim(dfClassClean)[2]
aiColOrder <- aiColOrder[!aiColOrder==iIndex]
aiColOrder <- c(iIndex,aiColOrder)

dfClassClean <- dfClassClean[,aiColOrder]
dfClassClean[,1] <- as.character(dfClassClean[,1])
astrMetaData <- c("Sample.ID","Gender","Country","Dataset")
for (strCat in astrMetaData){
  dfClassClean[strCat,1] <- strCat
}

dfClassClean["Sample.ID",1]
dfClassClean[,1] <- sub("_"," ",dfClassClean[,1],fixed=T )
###############################################################################################
# Write out the clean results
###############################################################################################

write.table(dfDetClean,opt$out_results,sep="\t",row.names=T,col.names=F,quote=F)
write.table(dfClassClean,opt$out_classresults,sep="\t",row.names=F,col.names=F,quote=F)
#write.table(matMeta,opt$out_meta,sep="\t",row.names=T,col.names=F)


