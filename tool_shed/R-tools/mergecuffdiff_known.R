# Merges sample fpkm with cuffdiff output file
# Maxime Caron - Dec 2011
# Mathieu Bourgey - Jan 2013
# args[1] = output dir
# args[2] = design file

args <- commandArgs(trailingOnly = TRUE)
output_dir<-args[1]
design_file<-args[2]

designs<-read.table(design_file, header=F, sep="\t", check.names=F)

# Iterate for each design

for(j in 2:ncol(designs)) {

	designName=designs[1,j] 
	reference<-read.table(paste(output_dir,"/known/",designName,"/isoform_exp.diff",sep=""), header=T, sep="\t")

	for(i in 2:nrow(designs)) {
        	sampleState=designs[i,j]
                if(sampleState !=0) {
                        sampleName=designs[i,1]
                        transcriptFile=paste(output_dir,"/known/",sampleName, "/transcripts.gtf",sep="")
			transcripts<-read.table(transcriptFile, header=T, sep="\t", quote='"')
			transcripts<-as.data.frame(transcripts[agrep("transcript", transcripts[,3]),9])
                        transcripts<-as.data.frame(matrix(unlist(strsplit(as.character(transcripts[,1]), ";")), nrow=nrow(transcripts), byrow=T))
                        ensemble_id<-as.data.frame(matrix(unlist(strsplit(as.character(transcripts[,2]), " ")), nrow=nrow(transcripts), byrow=T))
                        ensemble_id<-as.data.frame(ensemble_id[,3])
                        fpkm<-as.data.frame(matrix(unlist(strsplit(as.character(transcripts[,3]), " ")), nrow=nrow(transcripts), byrow=T))
                        finalSample<-cbind(ensemble_id, fpkm[,3])
			colnames(finalSample)=c(paste("id.",sampleName,sep=""),paste("fpkm.",sampleName,sep=""))
			reference<-merge(reference, finalSample, by.x=1, by.y=1)
		}
	}

write.table(reference, paste(args[1],"/known/", designs[1,j], "/isoform_exp.diff.with.fpkm.csv", sep=""), quote=F, row.names=F, sep="\t")

}
