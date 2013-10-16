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
## Effects by functional class
## 
## Count by effects
##  * deleting:
##     CODON_CHANGE_PLUS_CODON_DELETION
##     CODON_CHANGE_PLUS_CODON_INSERTION
##     CODON_DELETION
##     CODON_INSERTION
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
write.table(t(summaryTable),paste(outputBaseName,"Summary.table.tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)

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
pheatmap(t(as.matrix(chR.ord)),cluster_cols =F,cluster_rows =F,fontsize = 14,main="Change rate by sample and by chromosome")
dev.off()






## Effects by impact
## 
effectIlist=strsplit(gsub("%","",gsub(" ","",scan(listFiles[grep(fileExtensionRetained[4],listFiles)],sep="\n",what='character'))),",")
effectITable=NULL
for (i in 2:length(effectIlist)){
	effectITable=c(effectITable,effectIlist[[i]][1:2])
}
write.table(t(effectITable),paste(outputBaseName,"Effects.by.impact.tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)

## Effects by functional class
## 
effectIlist=strsplit(gsub("%","",gsub(" ","",scan(listFiles[grep(fileExtensionRetained[4],listFiles)],sep="\n",what='character'))),",")
effectITable=NULL
for (i in 2:length(effectIlist)){
	effectITable=c(effectITable,effectIlist[[i]][1:2])
}
write.table(t(effectITable),paste(outputBaseName,"Effects.by.impact.tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)

