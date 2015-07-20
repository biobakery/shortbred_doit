#! /usr/bin/env Rscript
# code obtained from 
# http://stackoverflow.com/questions/13258454/marking-specific-tiles-in-geom-tile-geom-raster
# and http://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/

# QUESTIONS
## Some classes are not listed in the gold standard. Do we assume those are 0? 


rm(list=ls())
library(ggplot2)
library(reshape)
library(optparse)



option_list <- list(
  make_option("--meta", default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\data\\ABR_Metadata\\FinalMetadataBrief.tab", type="character",
              help="Enter the metadata for the SB Markers."),
  make_option("--infants", default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\figures\\S4-InfantMetagenome\\WASHU_Contigs_mergedresults.tab", type="character",
              help="Infant Results."),
  make_option("--gs", default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\data\\InfantMetagenomes\\GoldStandard.txt", type="character",
              help="Gold Standard."),
  make_option("--out", default="C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-04-06\\InfantMetagenomes-Contigs.pdf", type="character",
              help="Gold Standard.")
)
opt <- parse_args(OptionParser(option_list=option_list))


#fileInfants <- "C:\\Users\\Jim\\Dropbox\\ShortBRED\\results\\2014-03-31\\InfantMetagenomes\\WASHU_mergedresults.tab"
fileInfants <- opt$infants
fileMetadata <- opt$meta
fileGS <- opt$gs


# Load the data, confirm that R knows it's numeric.
dfMeta <- read.table(fileMetadata,sep="\t",header=T)
colnames(dfMeta) <- c("ABR","Class")
dfInfants <- read.table(fileInfants,sep="\t",header=T,row.names=1)
print("Confirm data is numeric...")
sapply(dfInfants, is.numeric)
dfInfants$ABR <- rownames(dfInfants)
print ("Merge the data...")
dfInfants<- merge(dfInfants,dfMeta,by="ABR")


# Code for checking specifici ABR's of interest.
# dfInfants[dfInfants$Class=="rRNA_Methyltransferase",grep("F04_.*",colnames(dfInfants),perl=TRUE)]


dfInfants$ABR <- NULL



# Melt the data - 
meltInfants <- melt(dfInfants,id.vars="Class")

meltInfants$variable <- gsub("Forward_", "", meltInfants$variable)
meltInfants$variable <- gsub("Reverse_", "", meltInfants$variable)
meltInfants$variable <- gsub(".fastq..Count", "", meltInfants$variable)
meltInfants$variable <- gsub(".faa..Count", "", meltInfants$variable)
meltInfants$variable <- gsub("_.*", "", meltInfants$variable,perl=TRUE)
meltInfants <- aggregate(meltInfants$value,by=list(meltInfants$Class,meltInfants$variable),sum) 

colnames(meltInfants) <- c("ABR","Sample","ABR_Level")
#meltInfants$ABR_Level <- log(meltInfants$ABR_Level)

# Plot heatmap
print("Attempt heatmap...")
meltInfants$ABR_Level[meltInfants$ABR_Level<0] <- 0
g <- ggplot(meltInfants) +
  geom_tile( aes(x=Sample,y=ABR,fill=ABR_Level)) + 
  scale_fill_gradient( low = "white", high = "red", na.value="black", name = "" )

################################################################################
# Load gold standard, melt it, plot squares.
dfGS <- read.table(fileGS,sep="\t",header=T)
print("Check GS:")
print(dfGS[1,])


astrGS_Classes <- dfGS$ID
meltInfants <- meltInfants[meltInfants$ABR %in% astrGS_Classes ,]

meltGS <- melt(dfGS,id.vars="ID") 
colnames(meltGS) <- c("ABR","Sample","value")



meltGS$Sample <-as.integer(as.factor(meltGS$Sample))
meltGS$ABR <-as.integer(as.factor(meltGS$ABR))
meltGS<-meltGS[meltGS$value==1,]

y <- ggplot(meltInfants) +
  geom_tile( aes(x=Sample,y=ABR,fill=ABR_Level)) + 
  scale_fill_gradient( low = "white", high = "red", na.value="black", name = "" ) +
  geom_rect(data=meltGS, size=1, fill=NA, colour="black",
          aes(xmin=Sample - 0.5, xmax=Sample + 0.5, ymin=ABR - 0.5, ymax=ABR + 0.5)) +
  labs(title="ABR Abundance in Infant Metagenomes") +  
  theme(plot.title = element_text(size = rel(1.2)),axis.text.x = element_text(size=rel(0.4))) 
  

y

ggsave(file=opt$out, width=7, height=4) 
