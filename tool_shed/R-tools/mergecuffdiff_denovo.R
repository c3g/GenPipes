# Merges sample fpkm with cuffdiff denovo output file
# Maxime Caron - Jan 2012

# Arguments; output dir and design file

args <- commandArgs(trailingOnly = TRUE)
output_dir<-args[1]
input_dir<-args[2]
designFile<-args[3]
designs<-read.table(designFile, header=F, sep="\t", check.names=F)

# Iterate for each design

for(j in 2:ncol(designs)) {

	designName=designs[1,j]
	print(paste("design: ",designName, sep=""))
	resultFile<-paste(output_dir,"/denovo/",designName,"/isoform_exp.diff",sep="")
	mergedFile<-paste(output_dir,"/denovo/",designName,"/formated.merged.gtf", sep="")
	result<-read.table(resultFile, header=T, sep="\t", stringsAsFactors=F)
	merged<-read.table(mergedFile, header=F, sep="\t", quote='"', stringsAsFactors=F)
	tconsmerged<-as.data.frame(matrix(unlist(strsplit(as.character(merged[,1]), " ")), nrow=nrow(merged), byrow=T))
	oId<-as.data.frame(matrix(unlist(strsplit(as.character(merged[,2]), " ")), nrow=nrow(merged), byrow=T))
	ref_id<-as.data.frame(matrix(unlist(strsplit(as.character(merged[,3]), " ")), nrow=nrow(merged), byrow=T))
	tconsmerged<-as.data.frame(tconsmerged[,3])
	oId<-as.data.frame(oId[,3])
	ref_id<-as.data.frame(ref_id[,3])
	classcode<-merged[,4]
	mergedFinal<-cbind(tconsmerged, oId, ref_id, classcode)
	colnames(mergedFinal)=c("transcript_id", "oId", "nearest_ref", "classcode")

	# Merge with result file (isoform_exp.diff)

	writeIt<-merge(result, mergedFinal, by.x=1, by.y=1)
	writeIt<- writeIt[,c(15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17)]
	
	# Merge with sample FPKM

	for(i in 2:nrow(designs)) {
	#for(i in 2:2) {
		sampleState=designs[i,j]
		if(sampleState !=0) {
			sampleName=designs[i,1]
			transcriptFile=paste(input_dir,"/denovo/",sampleName, "/transcripts.gtf",sep="")
			transcripts<-read.table(transcriptFile, header=F, sep="\t", quote='"', stringsAsFactors=F)
			transcripts<-as.data.frame(transcripts[agrep("transcript", transcripts[,3]),9])
			transcripts<-as.data.frame(matrix(unlist(strsplit(as.character(transcripts[,1]), ";")), nrow=nrow(transcripts), byrow=T))
			tcons<-as.data.frame(matrix(unlist(strsplit(as.character(transcripts[,2]), " ")), nrow=nrow(transcripts), byrow=T))
			tcons<-as.data.frame(tcons[,3])
			fpkm<-as.data.frame(matrix(unlist(strsplit(as.character(transcripts[,3]), " ")), nrow=nrow(transcripts), byrow=T))
			finalSample<-cbind(tcons, fpkm[,3])
			colnames(finalSample)=c(paste("id.",sampleName,sep=""),paste("fpkm.",sampleName,sep=""))
			writeIt<-merge(writeIt, finalSample, by.x=1, by.y=1)
			#print(head(writeIt))
		}
	}
	write.table(writeIt, paste(args[1],"/denovo/", designs[1,j], "/isoform_exp.diff.with.fpkm.csv", sep=""), quote=F, row.names=F, sep="\t")
	
}
