#!/usr/bin/env Rscript

# Function that takes alpha diversity tables from qiime, and outdir and a prefix
# To write plot files.
makeplots <- function(mappingFile, rarefactionTables, outdir, prefix) {

  library(ggplot2)

  id    = NULL
  names = NULL
  df    = NULL
  xaxis = NULL

  fh2 = readLines(rarefactionTables)

  for(i in 1:length(fh2)){
    if(grepl("xaxis: ",fh2[i]) ){
      xaxis = gsub("xaxis: ", "", fh2[i])
      xaxis = as.numeric(unlist(strsplit(xaxis, "\t")))
    }
  
    if(grepl(">>",fh2[i]) ){
      id = gsub(">> ", "", fh2[i])
    }
    
    if(grepl("series",fh2[i]) ){
      row = gsub("series ", "", fh2[i])
      row = strsplit(row, "\t")
      value = as.numeric(row[[1]])
      df= rbind(df, data.frame(sample=id,value=value,seqEff=xaxis))
    }
  }

## Now read mapping file and merge the 2 data frames
#map = read.table(mappingFile, sep="\t",  comment.char = "#", header=F)
#map = data.frame(map)
#names(map) = c("sampleName", "sample")
#currDf = merge(map, df, by="sample")

p <- ggplot(df,aes(x=seqEff,y=value)) + 
  geom_point(aes(colour="#1E90FF")) + 
  facet_wrap(~sample) + 
  xlab("Sequencing effort") +
  ylab("Observed OTUs") +
  guides(colour=FALSE) +
  theme(
    panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
    axis.text.y=element_text(size=7, colour="black"),
    axis.text.x=element_text(size=10, colour="black", angle=90),
    axis.title=element_text(size=12),
    panel.background=element_blank()
  ) 
  pdf(file=paste0(outdir, "/", prefix, ".pdf"))
  print(p)
  dev.off()
  jpeg(file=paste0(outdir, "/", prefix, ".jpeg"), height=16, width=16, units="in", res=500)
  print(p)
  dev.off()
} 

usage=function(errM) {
        cat("\nUsage : Rscript pacBioAssemblyPlots.R [option] <Value>\n")
        cat("       -m        : mapping file\n")
        cat("       -r        : rarefaction file\n")
        cat("       -o        : outdir\n")
        cat("       -p        : prefix\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-m") {
		mappingFile=ARG[i+1]
	}else if (ARG[i] == "-r") {
		rarefactionTables=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}else if (ARG[i] == "-p") {
	  prefix=ARG[i+1]
	}
}

makeplots(mappingFile, rarefactionTables, outdir, prefix)

