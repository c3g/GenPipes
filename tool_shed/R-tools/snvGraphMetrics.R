## merge all sample readStat
## Mathieu bourgey 
## 2013/01/31

##load libraries
library(pheatmap)

args=commandArgs(TRUE)
listFiles=scan(args[1],sep="\n",what='character')
outputBaseName=args[2]

## what extension should be included in the listFiles
## Summary.table.csv
## Change.rate.by.chromosome.csv
## Changes.by.type.csv
## Effects.by.impact.csv
## Effects.by.functional.class.csv
## Count.by.effects.csv
## Count.by.genomic.region.csv
## Quality.csv
## Coverage.csv
## InDel.lengths.csv
## Base.changes.csv
## TsTv.summary.csv
## TsTv.All.variants.csv
## TsTv.Known.variants.csv
## Allele.frequency.csv
## Allele.frequency.All.variants.csv
## Allele.frequency.Known.variants.csv
## Codon.change.table.csv
## Amino.acid.change.table.csv
## Chromosome.change.table.csv
## changeRate.tsv
##

## Count by genomic region
## 
## Quality
##   * change it as a cumulative sum
## 
## Coverage
##   * change it as a cumulative sum and bin the coverage value
## 
## 
## InDel lengths
##   * barplot
## 
## Base changes
##   * 3d barplot
## 
## Ts/Tv : All variants
## Ts/Tv : Known variants
##   * pull the 2 rate together and plot them by samnple as barplot
## 
## Codon change table
##   * 3d barplot
## 
## Amino acid change table
##   * barplot
## 
## Chromosome change table
##   * supllementary figures in a zip

fileExtensionRetained=c("Summary.table.csv","TsTv.summary.csv","changeRate.tsv","Effects.by.impact.csv","Effects.by.functional.class.csv","Count.by.effects.csv","Count.by.genomic.region.csv","Quality.csv","Coverage.csv","InDel.lengths.csv","Base.changes.csv","TsTv.All.variants.csv","TsTv.Known.variants.csv","Codon.change.table.csv","Amino.acid.change.table.csv","Chromosome.change.table.csv")

## Summary table
##   * Number_of_variants_before_filter
##   * Number_of_variants_filtered_out
##   * Number_of_not_variants
##   * Number_of_variants_processed
##   * Number_of_known_variants
## Ts/Tv summary
##  * should be add to the summary table

sumT=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[1],listFiles)],sep="\n",what='character')),",")
valueSum=c("Number_of_variants_before_filter","Number_of_variants_filtered_out","Number_of_not_variants","Number_of_variants_processed","Number_of_known_variants")
summaryTable=NULL
for (i in 1:length(valueSum)){
	summaryTable=rbind(summaryTable,sumT[[grep(valueSum[i],sumT)]][1:2])
}
summaryTable=gsub("<br>(i.e.non-emptyID)","",summaryTable,fixed=T)
perS=list(c("%",as.character(round((as.numeric(summaryTable[2,2])/as.numeric(summaryTable[1,2]))*100,2))),c("%",as.character(round((as.numeric(summaryTable[3,2])/as.numeric(summaryTable[1,2]))*100,2))),c("%",as.character(round((as.numeric(summaryTable[5,2])/as.numeric(summaryTable[4,2]))*100,2))))
summaryTable=rbind(summaryTable[1:2,],perS[[1]],summaryTable[3,],perS[[2]],summaryTable[4:5,],perS[[3]])
sumTs=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[2],listFiles)],sep="\n",what='character')),",")
for ( i in 1:length(sumTs)) {
	summaryTable=rbind(summaryTable,sumTs[[i]][1:2])
}
colnames(summaryTable)=c("Summary_stats","Value")
write.table(t(summaryTable),paste(outputBaseName,"SummaryTable.tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)

## change rate
##   * graphs
##
chR=read.table(listFiles[grep(fileExtensionRetained[3],listFiles)],header=T,check.names=F,row.names=1)
cn=as.numeric(gsub("chr","",rownames(chR)))
cnpos=1:length(cn)
cnNu=order(cn[cnpos[!(is.na(cn))]])
cnNuV=chR[!(is.na(cn)),]
cnCh=order(rownames(chR)[cnpos[is.na(cn)]])
cnChV=chR[is.na(cn),]
chR.ord=rbind(cnNuV[cnNu,],cnChV[cnCh,])
jpeg(paste(outputBaseName,"changeRate.jpeg",sep="."),800,800)
pheatmap(t(as.matrix(chR.ord)),cluster_cols =F,cluster_rows =F,main="Change rate by sample and by chromosome",fontsize_row=6,fontsize_col=10)
dev.off()
pdf(paste(outputBaseName,"changeRate.pdf",sep="."),title="Change rate by sample and by chromosome",pointsize=5,paper='special')
pheatmap(t(as.matrix(chR.ord)),cluster_cols =F,cluster_rows =F,main="Change rate by sample and by chromosome",fontsize_row=3,fontsize_col=8)
dev.off()
write.table(t(as.matrix(chR.ord)),paste(outputBaseName,"changeRate.tsv",sep="."),quote=F,row.names=T,col.names=T,sep="\t")


## Effects by impact
## 
effectIlist=strsplit(gsub("%","",gsub(" ","",scan(listFiles[grep(fileExtensionRetained[4],listFiles)],sep="\n",what='character'))),",")
effectITable=NULL
nameEffect=NULL
for (i in 2:length(effectIlist)){
	effectITable=c(effectITable,effectIlist[[i]][2])
	nameEffect=c(nameEffect,effectIlist[[i]][1])
}
effectITable=rbind(nameEffect,effectITable)
write.table(effectITable,paste(outputBaseName,"EffectsImpact.tsv",sep="."),sep="\t",col.names=T,row.names=F,quote=F)

## Effects by functional class
## 
effectFlist=strsplit(gsub("%","",gsub(" ","",scan(listFiles[grep(fileExtensionRetained[5],listFiles)],sep="\n",what='character'))),",")
effectFTable=NULL
nameEffect=NULL
for (i in 2:length(effectFlist)){
	effectFTable=c(effectFTable,effectFlist[[i]][2])
	nameEffect=c(nameEffect,effectFlist[[i]][1])
}
effectFTable=rbind(nameEffect,effectFTable)
write.table(effectFTable,paste(outputBaseName,"EffectsFunctionalClass.csv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)

## Count by effects
##  * deleting:
##     CODON_CHANGE_PLUS_CODON_DELETION
##     CODON_CHANGE_PLUS_CODON_INSERTION
##     CODON_DELETION
##     CODON_INSERTION
## 

countByEffect=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[6],listFiles)],sep="\n",what='character')),",")
removeINDEL=c("CODON_CHANGE_PLUS_CODON_DELETION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_INSERTION")
countTable=NULL
countName=NULL
for (i in 2:length(countByEffect)){
	if (!(countByEffect[[i]][1] %in% removeINDEL)){
		countTable=c(countTable,countByEffect[[i]][2])
		countName=c(countName,countByEffect[[i]][1])
	}
}
countTableF=rbind(countName,countTable)
write.table(countTableF,paste(outputBaseName,"CountEffects.csv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
jpeg(paste(outputBaseName,"CountEffects.jpeg",sep="."),800,800)
par(las=2)
par(oma=c(10,4,4,1))
barplot(as.numeric(countTable),names.arg=countName,col=rainbow(length(countTable)),main="Total number of variant by effect type")
dev.off()
pdf(paste(outputBaseName,"CountEffects.pdf",sep="."),title="Total number of variant by effect type",paper='special')
par(las=2)
par(oma=c(10,4,4,1))
barplot(as.numeric(countTable),names.arg=countName,col=rainbow(length(countTable)),main="Total number of variant by effect type")
dev.off()




## Count by genomic region
## 

countByRegion=strsplit(gsub(" ","",scan(listFiles[grep(fileExtensionRetained[7],listFiles)],sep="\n",what='character')),",")
regionOrder=c("UPSTREAM","UTR_5_PRIME","SPLICE_SITE_ACCEPTOR","EXON","SPLICE_SITE_DONOR","INTRON","UTR_3_PRIME","DOWNSTREAM","INTERGENIC")
regionName=c("Up","5'","Splicing acceptor","Exon","Splicing donor","intron","3'","Down","Intergenic")
countTable=NULL
countName=NULL
for (i in 2:length(countByRegion)){
		countTable=c(countTable,countByRegion[[i]][2])
		countName=c(countName,countByRegion[[i]][1])
}
orderR=match(regionOrder,countName)
countTableF=rbind(countName[orderR],countTable[orderR])
write.table(countTableF,paste(outputBaseName,"CountRegions.csv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
jpeg(paste(outputBaseName,"CountRegions.jpeg",sep="."),800,800)
par(las=2)
par(oma=c(10,4,4,1))
barplot(as.numeric(countTable),names.arg=regionName,col=rainbow(length(countTable)),main="Total number of variant by region type")
dev.off()
pdf(paste(outputBaseName,"CountRegions.pdf",sep="."),title="Total number of variant by region type",paper='special')
par(las=2)
par(oma=c(10,4,4,1))
barplot(as.numeric(countTable),names.arg=regionName,col=rainbow(length(countTable)),main="Total number of variant by region type")
dev.off()
