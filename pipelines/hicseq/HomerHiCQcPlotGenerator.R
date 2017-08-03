############################################################################
##                                                                        ##
##                              Rola Dali                                 ##  
##                        Homer HiC QC autoPlot                           ##
##                              April 2016                                ##
##                                                                        ##
############################################################################

## Rscript HomerHiCQcPlotGenerator.R <SampleName> <WorkingDirectory> <OutputDirectory>

## Script that generates QC plots from Homer HiC pipeline:

#petagLocalDistribution.txt: shows the relationship between the 5' ends of the paired reads; used to determine fragment size of the Hi-C fragments.

#petagDistDistribution.txt: histogram with fraction of paired-end reads that are found at different distances from one another.  The file stops at 300 million bp, and is divided into 1kb bins.  The top of the file also contains the fraction of interchromosomal interactions.

#petagRestrictionDistribution.<Res>.mis<0>.txt: shows the distribution of reads around the restriction site. mis<0> shows allowed mismatched in case of star activity.

#tagCountDistribution.txt: File contains a histogram of clonal read depth, showing the number of reads per unique position. --> duplication level

#tagAutocorrelation.txt: autocorrelation routine creates a distribution of distances between adjacent reads in the genome. similar to petagLocalDistribution?

#tagInfo.txt: number of tags per chr

####################################################################################

## Rscript HomerHiCQcPlotGenerator.R <SampleName> <WorkingDirectory> <OutputDirectory>



## TODO: add plot for petagDistDistribution.txt


## command line arguments for R script:
args<-commandArgs(TRUE)
sampleName <- args[1]
workingDir <- args[2]

if(!is.na(args[3])){
  outputDir <- args[3]
} else {
  outputDir <- "HomerQcPlots"
}

print(workingDir)

setwd(workingDir)
if (!file.exists(outputDir)) {dir.create(outputDir)}

print("... ... Loading libraries")

## load libraries:
library(ggplot2)
library(reshape)
library(stringr)

## read input files and plot:

print("... ... producing petagLocalDistribution.pdf")

#petagLocalDistribution.txt
file <- list.files(pattern = "petag.LocalDistribution.txt", full = TRUE)

if (length(file) != 0 ){
if(length(file) >1){
  file <- list.files(workingDir, pattern = paste(sampleName ,"_petag.LocalDistribution.txt", sep=""), full = TRUE)
}

data <- read.table(file, header=T, sep="\t")

data <- melt(data = data, id="Local.Distance.in.bp.between.PE.tags")
colnames(data) <- c("Local.Distance.in.bp.between.PE.tags", "strand", "value")

pdf(file = paste(outputDir, "/", sampleName, "_petag.LocalDistribution.pdf", sep=""), width = 14)
print(ggplot() + geom_line(data=data, aes(x=Local.Distance.in.bp.between.PE.tags, y=value , color=strand), size=1.5)  + xlim (-1000,1000) + xlab("Distance from 5end of First Read") + ylab("Counts per bp") + ggtitle(paste(sampleName, "_petag.LocalDistribution", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())))
dev.off() } else { print ("file not found")}

print("... ... producing tagLengthDistribution.pdf")

#tagLengthDistribution.txt
file <- list.files(pattern = "tagLengthDistribution.txt", full = TRUE)
#print(file)
if (length(file) != 0 ){
if(length(file) >1){
  file <- list.files(workingDir, pattern = paste(sampleName ,"_tagLengthDistribution.txt", sep=""), full = TRUE)
}

data <- read.table(file, header=T, sep="\t")
colnames(data) <- c("TagLength",  "FractionOfTags")

pdf(file = paste(outputDir, "/",sampleName, "_tagLengthDistribution.pdf", sep=""), width = 14, height = 7)
print(ggplot() + geom_line(data=data, aes(x=TagLength, y=FractionOfTags), color="red", size=1.5) + xlab(" Tag length") + ylab("Fraction of tags") + ggtitle(paste(sampleName, "_tagLengthDistribution", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())))
dev.off() } else {print("file not found")}

print("... ... producing petagRestrictionDistribution.pdf")


#petagRestrictionDistribution.<Res>.mis<0>.txt
file <- list.files(pattern = "petagRestrictionDistribution", full = TRUE)

if (length(file) != 0 ){
if(length(file) >1){
  file <- list.files(workingDir, pattern = paste(sampleName ,"_petagRestrictionDistribution", sep=""), full = TRUE)
}

data <- read.table(file, header=T, sep="\t")
colnames(data) <- c("DistanceFromResSite", "+Strand", "-Strand")
data <- melt(data, id="DistanceFromResSite")
colnames(data) <- c("DistanceFromResSite", "strand", "value")

cutSite <- str_match(string = file, pattern = "petagRestrictionDistribution.*..mis")
cutSite <- gsub(x = gsub(x= cutSite, pattern = "petagRestrictionDistribution.", replacement = ""), pattern = ".mis", replacement = "")

pdf(file = paste(outputDir, "/", sampleName, "_petagRestrictionDistribution.pdf", sep=""), width = 14, height = 7)
print(ggplot() + geom_line(data=data, aes(x=DistanceFromResSite, y=value, color=strand), size=1.5) + xlab("Distance From Restriction Site - ") + ylab("Fraction of tags") + ggtitle(paste(sampleName, "_tagLengthDistribution_", cutSite, sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())))
dev.off() } else {print("file not found")}

print("... ... producing tagCountDistribution.pdf")

#tagCountDistribution.txt

file <- list.files(pattern = "tagCountDistribution.txt", full = TRUE)

if (length(file) != 0 ){
if(length(file) >1){
  file <- list.files(workingDir, pattern = paste(sampleName ,"_tagCountDistribution.txt", sep=""), full = TRUE)
}

data <- read.table(file, header=T, sep="\t")
colnames(data) <- c("TagsPerPosition", "FractionOfTags")

pdf(file = paste(outputDir, "/", sampleName, "_tagCountDistribution.pdf", sep=""), width = 14, height = 7)
print(ggplot(data=data, aes(TagsPerPosition, FractionOfTags)) + geom_bar(stat = "identity", fill="grey")  + xlab("Reads per position") + ylab("Fraction of total reads") + ggtitle(paste(sampleName, "_tagCountDistribution", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())) +  scale_x_continuous(limits=c(0.5, 10), breaks=1:10) )
dev.off() } else {print("file not found")}

print("... ... producing tagAutocorrelation.pdf")

#tagAutocorrelation.txt


file <- list.files(pattern = "tagAutocorrelation.txt", full = TRUE)

if (length(file) != 0 ){
if(length(file) >1){
  file <- list.files(workingDir, pattern = paste(sampleName ,"_tagAutocorrelation.txt", sep=""), full = TRUE)
}

data <- read.table(file, header=T, sep="\t")
colnames(data) <- c("DistanceInBp", "+Strand", "-Strand")
data <- melt(data, id="DistanceInBp")
colnames(data) <- c("DistanceInBp", "strand", "value")


pdf(file = paste(outputDir, "/",sampleName, "_tagAutocorrelation.pdf", sep=""), width = 14, height = 7)
print(ggplot() + geom_line(data=data, aes(x=DistanceInBp, y=value, color=strand), size=1.5) + xlab("Relative distance btn Reads (bp)") + ylab("Total read pairs") + ggtitle(paste(sampleName, "_tagAutocorrelation", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())))
dev.off() } else { print("file not found")}

print("... ... producing tagInfo.pdf")

#tagInfo.txt 
file <- list.files(pattern = "tagInfo.txt", full = TRUE)

if (length(file) != 0 ){
if(length(file) >1){
  file <- list.files(workingDir, pattern = paste(sampleName ,"_tagInfo.txt", sep=""), full = TRUE)
}

data <- read.table(file, header=T, sep="\t", skip = 11)
colnames(data) <- c("chr", "UniquePositions", "TotalTags")
chrList <- paste("chr" , c( seq(1:22), "X", "Y"),sep="" )
data <- data[data$chr %in% chrList,]
data = melt(data, id="chr")

pdf(file = paste(outputDir, "/",sampleName, "_tagInfo.pdf", sep=""), width = 14, height = 7)
print(ggplot() + geom_point(data=data, aes(x=chr, y=value, color=variable)) + xlab("chr") + ylab("Tag Counts") + ggtitle(paste(sampleName, "_tagAutocorrelation", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())) + scale_x_discrete(limits=chrList) )
dev.off() }  else { print("file not found")}

## if GCcheck is Run, 3 more files would be produced. Check and plot:


#tagFreq.txt 
file <- list.files(pattern = "tagFreq.txt", full = TRUE)

if (length(file) != 0 ){
print("... ... producing tagFreq.txt")

if(length(file) >1){
  file <- list.files(pattern = paste(sampleName ,"_tagFreq.txt", sep=""), full = TRUE)
}

data <- read.table(file, header=T, sep="\t")
data <- data[,1:5]
data <- melt(data, id="Offset")
colnames(data) <- c("bp", "nucleotide", "value")

pdf(file = paste(outputDir, "/",sampleName, "_tagFreq.pdf", sep=""), width = 14, height = 7)
print(ggplot() + geom_line(data=data, aes(x=bp, y=value, color=nucleotide)) + xlab("Distance from 5end of Read") + ylab("Nucleotide Frequency") + ggtitle(paste(sampleName, "_tagFreq", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())) )
dev.off()
}



#tagGCcontent.txt 
file <- list.files(pattern = "tagGCcontent.txt", full = TRUE)
genomeFile <- list.files(pattern = "genomeGCcontent.txt", full = TRUE)

if (length(file) != 0 ){
  print("... ... producing tagGCcontent.pdf")
  
  if(length(file) >1){
    file <- list.files(workingDir, pattern = paste(sampleName ,"_tagGCcontent.txt", sep=""), full = TRUE)
    genomeFile <- list.files(workingDir, pattern = paste(sampleName ,"_genomeGCcontent.txt", sep=""), full = TRUE)
  }
  
  data <- read.table(file, header=T, sep="\t")
  colnames(data) <- c("GC", "Total", "Normalized")
  data$sample <- "sample"
  
  genome <- read.table(genomeFile, header=T, sep="\t")
  colnames(genome) <- c("GC", "Total", "Normalized")
  genome$sample <- "genome"
  
  data <- rbind(data, genome)
  
  
  pdf(file = paste(outputDir, "/",sampleName, "_tagGCcontent.pdf", sep=""), width = 14, height = 7)
  print(ggplot() + geom_line(data=data, aes(x=GC, y=Normalized, color = sample), size=1.5) + xlab("GC content of fragments(%)") + ylab("Normalized Fraction") + ggtitle(paste(sampleName, "_tagGCcontent", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())) + scale_x_continuous(breaks=c(seq(from=0, to=1, by = 0.1)))) 
  dev.off()
}

print("... ... Script Complete")


