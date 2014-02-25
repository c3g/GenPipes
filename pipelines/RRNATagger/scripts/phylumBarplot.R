#!/usr/bin/env Rscript

# Function that takes a infile: (filtered_summary.csv) and an outdir where to
# To write plot files.
generatePlots <- function(infile, outdir, prefix) {
	library(ggplot2)
	library(scales)
	library(reshape2)
	library(grid)
  
	print(paste0("[DEBUG] infile: ",infile))
	print(paste0("[DEBUG] outdir: ",outdir))
	print(paste0("[DEBUG] prefix: ",prefix))
	
	vColors = c(
	  '#00FF00',
	  '#FF8080',
	  '#FF00FF',
	  '#0000FF',
	  '#00CCFF',
	  '#CCFFFF',
	  '#CCFFCC',
	  '#99CCFF',
	  '#CC99FF',
	  '#FFCC99',
	  '#3366FF',
	  '#33CCCC',
	  '#99CC00',
	  '#FF99CC',
	  '#FFCC00'
	)

	tData <- read.csv(file=infile, header=T, sep="\t")
	tData2 <- cbind(tData, rowSums(tData[2:ncol(tData)]))
	tData3 <- tData2[order(-tData2[, ncol(tData2)]),]
	tData3[, ncol(tData3)] <- NULL
	tData4 <- tData3[1:15,]
	tData5 <- melt(tData4)
 
	print("[DEBUG] Loaded data...")

	outfile1 <- paste0(outdir, prefix, ".pdf")
	outfile2 <- paste0(outdir, prefix, ".jpeg")
 
	p <- ggplot(data=tData5, aes(order=order(-variable), x = variable, y=value, fill=Taxon)) + 
	  geom_bar(stat="identity") + 
	    theme(
	      panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
	      axis.text.x=element_text(size=7, colour="black", angle=90), 
	      axis.title=element_text(size=16),
	      panel.grid.major=element_line(colour="black", linetype="dotted"),
	      panel.grid.minor=element_blank(),
	      panel.background=element_blank(),
	 	  legend.key.size = unit(0.45, "cm"),
		  legend.margin = unit(1, "cm")
	    ) +  scale_fill_manual(values=vColors)
	pdf( file=outfile1, height=8, width=16 )
	print(p)
	dev.off()

	jpeg( file=outfile2, height=8, width=16, units="in", res=500)
	print(p)
	dev.off()
	
	print("[DEBUG] Printed plots...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript phylumBarplot.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outdir\n")
        cat("       -p        : prefix\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if (ARG[i] == "-i") {
		infile=ARG[i+1]
	} else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}else if (ARG[i] == "-p") {
		prefix=ARG[i+1]
	}
}

generatePlots(infile, outdir, prefix)

