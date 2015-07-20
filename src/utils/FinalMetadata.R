#! /usr/bin/env Rscript
rm(list=ls()) 
library(optparse)

file.sep <- .Platform$file.sep
dirARDB <- paste("C:","Users","Jim","Dropbox","ShortBRED","data","ARDB",sep=file.sep)
dirWUSTL <- paste("C:","Users","Jim","Dropbox","ShortBRED","data","WashU_ABR",sep=file.sep)
dirWD <- paste("C:","Users","Jim","Dropbox","ShortBRED","data","ABR_Metadata",sep=file.sep) 

option_list <- list(
  make_option("--ardb", default=paste(dirARDB,"ARDBmeta.txt",sep=file.sep), type="character",
              help="ARDB metadata file."),
  make_option("--combined", default=paste(dirWD,"Metadata-DetailedWUSTL.txt",sep=file.sep), type="character",
              help="Combined metadata file"),
  make_option("--fams", default=paste(dirWD,"SB_Fams.tab",sep=file.sep), type="character",
              help="final.map file produced by shortbred-identify"),
  make_option("--out", default=paste(dirWD,"FinalMetadata.tab",sep=file.sep), type="character",
              help="Combined metadata file with markers"),
  make_option("--small", default=paste(dirWD,"FinalMetadataBrief.tab",sep=file.sep), type="character",
              help="Combined metadata file with markers")
)

opt <- parse_args(OptionParser(option_list=option_list))

##################################################################################
# Load Data
##################################################################################


dfARDB <- read.table(opt$ardb, header=TRUE, sep=";", row.names=1,
                     as.is=TRUE)
dfMeta <- read.table(opt$combined, header=TRUE, sep="\t",
                      as.is=TRUE,quote="")
dfMap <- read.table(opt$fams, header=FALSE, sep="\t",
                     as.is=TRUE)



##################################################################################
# Format WUSTL
# Take the final.map file, 
colnames(dfMap)<- c("Family","Marker")
dfMap$Marker <- NULL


dfFams <- as.data.frame(unique(dfMap$Family))
colnames(dfFams) <- "Family"

dfFams$WUSTL <- gsub("_[0-9]*$", "", dfFams$Family, perl =TRUE)

astrKeep <- c("ID_Orig","Class","Description","Source")

#################################################################################
# merge on ARDB 
rownames(dfFams) <- dfFams$Family
dfCombined <- merge(dfFams,dfARDB,by="row.names",all.x=TRUE)
dfCombined <- dfCombined[,2:4]
colnames(dfCombined)[3] <- "Class"


################################################################################
# replace NA's with WUSTL

dfCombined$Class[is.na(dfCombined$Class)] <- dfCombined[is.na(dfCombined$Class),"WUSTL"]

#################################################################################

dfCombined <- dfCombined[,c(1,3)]

colnames(dfCombined) <- c("Family","ID_Orig")
dfCombined[dfCombined$Family=="NP_752857",]
dfMeta[dfMeta$ID_Orig=="mdfa",]

###################################################################################
# Prep dfMeta


dfFinal <- merge(dfCombined,dfMeta,by="ID_Orig")
dfFinal <- rbind(dfFinal,c("unassigned","unassigned","unassigned"))

dfFinal[dfFinal$Family=="NP_752857",]
dfFinal <- dfFinal[!is.na(dfFinal$Family) & !is.na(dfFinal$Merged.ID),]

dfFinal$Merged.ID <- gsub("_", " ", dfFinal$Merged.ID, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

write.table(dfFinal,sep="\t",file=opt$out,row.names=FALSE,fileEncoding="UTF-8")

dfFinal <- dfFinal[,c("Family","Merged.ID")]
write.table(dfFinal,sep="\t",file=opt$small,row.names=FALSE,fileEncoding="UTF-8")
