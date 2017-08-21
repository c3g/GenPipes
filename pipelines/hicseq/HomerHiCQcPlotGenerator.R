
################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################


## Rscript HomerHiCQcPlotGenerator.R <SampleName> <WorkingDirectory> <OutputDirectory>

## Script that generates QC plots from Homer HiC pipeline:

#petagLocalDistribution.txt: shows the relationship between the 5' ends of the paired reads; used to determine fragment size of the Hi-C fragments.

#petagDistDistribution.txt: histogram with fraction of paired-end reads that are found at different distances from one another.  The file stops at 300 million bp, and is divided into 1kb bins.  The top of the file also contains the fraction of interchromosomal interactions.

#petagRestrictionDistribution.<Res>.mis<0>.txt: shows the distribution of reads around the restriction site. mis<0> shows allowed mismatched in case of star activity.

#tagCountDistribution.txt: File contains a histogram of clonal read depth, showing the number of reads per unique position. --> duplication level

#tagAutocorrelation.txt: autocorrelation routine creates a distribution of distances between adjacent reads in the genome. similar to petagLocalDistribution?

#tagInfo.txt: number of tags per chr

#petagDistDistribution.txt histogram describing the fraction of paired-end reads that are found at different distances from one another.  The file stops at 300 million bp, and is divided into 1kb bins.

####################################################################################

## Rscript HomerHiCQcPlotGenerator.R <SampleName> <WorkingDirectory> <OutputDirectory>


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

  data <- read.table(file, header=T, sep="\t")
  
  data <- melt(data = data, id="Local.Distance.in.bp.between.PE.tags")
  colnames(data) <- c("Local.Distance.in.bp.between.PE.tags", "strand", "value")
  
  pdf(file = paste(outputDir, "/", sampleName, "_petag.LocalDistribution.pdf", sep=""), width = 14)
  print(ggplot() + geom_line(data=data, aes(x=Local.Distance.in.bp.between.PE.tags, y=value , color=strand), size=1.5)  + xlim (-1000,1000) + xlab("Distance from 5end of First Read") + ylab("Counts per bp") + ggtitle(paste(sampleName, "_petag.LocalDistribution", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())))
  dev.off() } else { print ("file not found")}

print("... ... producing tagLengthDistribution.pdf")

#tagLengthDistribution.txt
file <- list.files(pattern = "tagLengthDistribution.txt", full = TRUE)

if (length(file) != 0 ){

  data <- read.table(file, header=T, sep="\t")
  colnames(data) <- c("TagLength",  "FractionOfTags")
  
  pdf(file = paste(outputDir, "/",sampleName, "_tagLengthDistribution.pdf", sep=""), width = 14, height = 7)
  print(ggplot() + geom_line(data=data, aes(x=TagLength, y=FractionOfTags), color="red", size=1.5) + xlab(" Tag length") + ylab("Fraction of tags") + ggtitle(paste(sampleName, "_tagLengthDistribution", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())))
  dev.off() } else {print("file not found")}

print("... ... producing petagRestrictionDistribution.pdf")


#petagRestrictionDistribution.<Res>.mis<0>.txt
file <- list.files(pattern = "petagRestrictionDistribution", full = TRUE)
  
if (length(file) != 0 ){
  
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
  
  data <- read.table(file, header=T, sep="\t")
  colnames(data) <- c("TagsPerPosition", "FractionOfTags")
  
  pdf(file = paste(outputDir, "/", sampleName, "_tagCountDistribution.pdf", sep=""), width = 14, height = 7)
  print(ggplot(data=data, aes(TagsPerPosition, FractionOfTags)) + geom_bar(stat = "identity", fill="grey")  + xlab("Reads per position") + ylab("Fraction of total reads") + ggtitle(paste(sampleName, "_tagCountDistribution", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())) +  scale_x_continuous(limits=c(0.5, 10), breaks=1:10) )
  dev.off() } else {print("file not found")}

print("... ... producing tagAutocorrelation.pdf")

#tagAutocorrelation.txt


file <- list.files(pattern = "tagAutocorrelation.txt", full = TRUE)

if (length(file) != 0 ){
  
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



print("... ... producing petagDistDistribution.pdf")

#petag.FreqDistribution_1000.txt

file <- list.files(pattern = "petag.FreqDistribution_1000.txt", full = TRUE)

if (length(file) != 0 ){
  
  data <- read.table(file, header=F, sep="\t")
  colnames(data) <- c("bins", "FractionOfTags")
  
  interChr <- data[1,2]
  interChr <- as.character(interChr)
  interChr <- gsub(pattern = "Fraction of total PE tags", replacement = "", x = interChr)
  
  data <- data[-1,]
  data[,1] <- as.character(data[,1])
  data[nrow(data),1] <- 300000000
  data[,1] <- as.numeric(data[,1])
  data[,2] <- as.numeric(as.character(data[,2]))
  
  png(file = paste(outputDir, "/", sampleName, "_petag.FreqDistribution_1000.png", sep=""), width = 14 , height = 7, units = "in", res=150)
  print(ggplot(data=data, aes(x= log10(bins), y=log10(FractionOfTags))) + xlab("Distance between Regions") + geom_point() + geom_line () + ylab( "Fraction of total PE tags") + ggtitle(paste(sampleName, "_petag.FreqDistribution_1000", sep="")) + theme_set(theme_bw(12) + theme_bw()+ theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))  + annotate("text", x=log10(data[nrow(data)/2,1]), y = log10(1000), label= interChr))
  dev.off() } else {print("file not found")}




print("... ... Script Complete")
