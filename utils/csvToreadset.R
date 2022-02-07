
################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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

## This script takes the csv project file downloadable from nanuq and creates the readset output tsv file needed to run the mugqic pipelines.

usage <- function (){
  print("Rscript csvToreadset.R csv_input.csv output_name [fastq|bam] data_path")
  print("   csv_input.csv is a comma separated file that contains experimental details about each sample. It can be downloaded from nanuq's Readset page (not Libraries page).")
  print("    indicate the name of the output readset file.")
  print("    indicate the type of input as fastq files or a bam.")
  print("    data_path is the relative path to the data from the project folder")
  
}


## parse command line arguments:

args<-commandArgs(TRUE)

if (length(args) != 4){
  usage()
  stop ("Error! missing arguments")
  
}

file <- args[1]
output <- args[2]
type <- args[3]
path <- args[4]


## open file and pull needed fields:
## needed fields are:
fields <- c("Sample", "Readset", "Library", "RunType", "Run", "Lane", "Adapter1", "Adapter2", "QualityOffset", "BED", "FASTQ1",  "FASTQ2",  "BAM")

## read file and filter columns:

data <- read.csv(file, as.is=T)

keep <- c("Name", "Library.Name", "Library.Barcode", "Run.Type", "Run", "Region", "Adaptor.Read.1..NOTE..Usage.is.bound.by.Illumina.Disclaimer.found.on.Nanuq.Project.Page.", "Adaptor.Read.2..NOTE..Usage.is.bound.by.Illumina.Disclaimer.found.on.Nanuq.Project.Page.", "Quality.Offset", "BED.Files", "Filename.Prefix")

keepFinal <- c(keep[-length(keep)], "1",  "2",  "3")

data <- data[, colnames(data) %in% keep]



fill_files <- function (data) {
  if (type == "fastq"){
  if (data[4] == "PAIRED_END"){
    fq1 <- file.path(path, paste0(data[11], "_R1.fastq.gz"))
    fq2 <- file.path(path, paste0(data[11], "_R2.fastq.gz"))
    bam <- ""
    
  } else if (data[4] == "SINGLE_END"){
    fq1 <-  file.path(path, paste0(data[11], "_R1.fastq.gz"))
    fq2 <- ""
    bam <- ""
    
  } else {
    stop ("Error! 'Run.Type' must be 'SINGLE_END' or 'PAIRED_END'")
  }
  
} else if (type == "bam") {
  fq1 <- ""
  fq2 <- ""
  bam <- file.path(path, paste0(data[11], ".bam"))
} else {
  stop ("Error! 'Type' must be 'bam' or 'fastq'")
  
}
  return (c(fq1,fq2,bam))
}


files <- t(apply(data, 1, function(x) fill_files(x)))
idx <- which(colnames(data) == "Filename.Prefix")
data <- data[,-idx]
data <- cbind(data, files)

data[is.na(data)] <- ""

data <- data[, keepFinal]

colnames(data) <- fields

write.table(data, output, sep="\t", row.names=F, quote = F)
