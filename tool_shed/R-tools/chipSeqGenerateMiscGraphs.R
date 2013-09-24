# Generates graphs based on different statistics
# Maxime Caron - Jan 2012

library(plotrix)

args <- commandArgs(TRUE)
design.file<-args[1]
output_dir<-args[2]
designs<-read.table(design.file, header=F, sep="\t", check.names=F)



for(i in 2:ncol(designs[1,])) {

	design<-unlist(strsplit(as.character(designs[1,i]),","))
	if(design[2] == "N") {
		designName<-design[1]

		# TSS categories stats
		postscript(paste("graphs/", designName, "_Misc_Graphs.ps", sep=""))
		par(mfrow=c(2,2))
		tss<-paste(output_dir, "/annotation/", designName, "/", designName, ".tss.stats.csv",sep="")
		d1<-read.table(tss, header=T, sep=",", check.names=F)
		for(i in 1:nrow(d1)) {
			slices <- c(d1[,1] + d1[,2], d1[,3], d1[,4], d1[,5], d1[,6], d1[,7])
			lbls <- c("gene", names(d1[3:length(d1)]))
			pct <- round(slices/sum(slices)*100)
			lbls <- paste(lbls, "(",pct, sep="") # add percents to labels
			lbls <- paste(lbls,"%)",sep="") # ad % to labels 
			pie(slices,labels=lbls, main=paste("Location analysis of binding sites\ndesign ", designName, sep=""))
		}
       		
		# Exon intron stats
		exons<-paste(output_dir, "/annotation/", designName, "/", designName, ".exon.stats.csv",sep="")
		if(file.info(exons)$size == 0) {
			d1=data.frame(c(0))
		}
		else {
			d1<-read.table(exons, header=F, sep=",", check.names=F)
		}
		introns<-paste(output_dir, "/annotation/", designName, "/", designName, ".intron.stats.csv",sep="")

		if(file.info(introns)$size == 0) {
                        d1=data.frame(c(0))
                }
                else {
                        d2<-read.table(introns, header=F, sep=",", check.names=F)
                }	
		hist(d1[,1], breaks=length(levels(as.factor(d1[,1]))), xlim=c(0,20), main=paste("Distribution of peaks found within exons\ndesign", designName, sep=""), xlab="Exon", ylab="Number of peaks")
		hist(d2[,1], breaks=length(levels(as.factor(d2[,1]))), xlim=c(0,20), main=paste("Distribution of peaks found within introns\ndesign", designName, sep=""), xlab="Intron", ylab="Number of peaks")
		
		# Distance to TSS

		distance<-paste(output_dir, "/annotation/", designName, "/", designName, ".tss.distance.csv",sep="")
        	d1<-read.table(distance, header=F, sep=",", check.names=F)
       		d1<-subset(d1, d1[,1] > -10000 & d1[,1] < 10000)
		hist(d1[,1], breaks=seq(-10000,10000,1000), main=paste("Distribution of peak distances relative to TSS\ndesign ", designName, sep=""), xlab="Distance to TSS (bp)", ylab="Number of peaks")
		
		dev.off()
	}

}

# Generate table stats with number of peaks, etc

toPrint<-rbind(c("design", "number of peaks", "percent near tss", "median peak height", "highest peak", "lowest peak", "avg peak width"))

for(i in 2:ncol(designs[1,])) {
        design<-unlist(strsplit(as.character(designs[1,i]),","))
        if(design[2] == "N") {
                designName<-design[1]
		file<-paste(output_dir, "/MACS/", designName, "_peaks.bed",sep="")
  		d1<-read.table(file, header=F, sep="\t", check.names=F)
       		nPeaks<-nrow(d1)
		averagePeakWidth<-round(mean(d1[,3]-d1[,2]),0)
		
		file<-paste(output_dir, "/MACS/", designName, "_summits.bed",sep="")
                d1<-read.table(file, header=F, sep="\t", check.names=F)
		highestPeak<-max(d1[,5])
		lowestPeak<-min(d1[,5])
		medianPeakHeight<-median(d1[,5])
		file<-paste(output_dir, "/annotation/", designName, "/", designName, ".tss.distance.csv",sep="")
                d1<-read.table(file, header=F, sep=",", check.names=F)
		d2<-subset(d1, d1[,1]>-1000 & d1[,1]<1000)
		percentNearTSS<-round(nrow(d2)/nrow(d1), 2)*100
		toPrint<-rbind(toPrint, c(designName, nPeaks, percentNearTSS, medianPeakHeight, highestPeak, lowestPeak, averagePeakWidth))
		
	}
}
#write.table(rbind(dName, nPeaks, percentNearTSS, medianPeakHeight, highestPeak, lowestPeak, averagePeakWidth), paste(output_dir, "/annotation/peak_stats.csv", sep=""), row.names=F, col.names=F, quote=F, sep=",")
write.table(toPrint, paste(output_dir, "/annotation/peak_stats.csv", sep=""), row.names=F, col.names=F, quote=F, sep=",")

