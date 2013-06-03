# Performs differential gene expression with DESeq
# Written by Maxime Caron - Dec 2011
# Usage : Rscript deseq.R -d path_design -c path_rawcountfile -o output_dir

library(DESeq)

# Usage

usage=function(errM) {
	cat("\nUsage : Rscript deseq.R [option] <Value>\n")
	cat("       -d        : design file\n")
	cat("       -c        : raw count file\n")
	cat("       -o        : output directory\n")
	cat("       -h        : this help\n\n")
	stop(errM)
}
set.seed(123456789)
perform_dge=function(counts, groups, count_limit, path) {

# Retain row which have > count_limit

counts<-counts[rowSums(counts) > count_limit,]

# Normalize and do test

cds<-newCountDataSet(counts, groups)
cds<-estimateSizeFactors(cds)
sizeFactors(cds)
if(length(groups)==2) {
cds<-estimateDispersions(cds, method="blind", sharingMode="fit-only")
}
else {
cds<-estimateDispersions(cds, method="pooled")
}

res<-nbinomTest(cds, "1", "2" )
res[,c(5,6)] = round(res[,c(5,6)], digits=3)
res[,7] = as.numeric(format(res[,7], digits=2))
res[,8] = as.numeric(format(res[,8], digits=2))
colnames(res)[c(1,7,8)] = c("id", "deseq.p-value", "deseq.adj.pvalue")
write.table(res[order(res[,8]), c(1,7,8)], paste(path,"deseq_results.csv",sep="/"), quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
fileOpen=paste(path,"edger_results.csv",sep="/")
d1<-read.table(fileOpen, header=T, sep="\t", quote="")
d2<-merge(d1, res[, c(1,7,8)], by.x=1, by.y=1, sep="\t")
d2<-d2[order(d2[,(ncol(d2)-1)]),]
vecWrite<-c(1:4, (ncol(d2)-1), ncol(d2), 5:6, 7:(ncol(d2)-2))
write.table(d2[,vecWrite], paste(path,"dge_results.csv",sep="/"), quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
}


##################################

ARG = commandArgs(trailingOnly = T)
## default arg values
count_limit=9
fpath="."
design_file=""
rawcount_file=""
out_path=""
## get arg variables
for (i in 1:length(ARG)) {
	if (ARG[i] == "-d") {
		design_file=ARG[i+1]
	} else if (ARG[i] == "-c") {
		rawcount_file=ARG[i+1]
	} else if (ARG[i] == "-o") {
		out_path=ARG[i+1]
	} else if (ARG[i] == "-h") {
		usage("")
	}
}
## check arg consitency
if (!(file.exists(design_file))) {
	usage("Error : Design file not found")
}
if (!(file.exists(rawcount_file))) {
	usage("Error : Raw count file not found")
}
if (out_path == "") {
	usage("Error : Output directory not specified")
}
tmpFP=strsplit(fpath,"")
if (tmpFP[[1]][length(tmpFP[[1]])] == "/" ) {
	bckS=""
} else {
	bckS="/"
}
tmpOP=strsplit(out_path,"")
if (tmpOP[[1]][length(tmpOP[[1]])] == "/" ) {
	out_path=paste(tmpOP[[1]][1:(length(tmpOP[[1]]-1))],collapse="")
}

design = read.csv2(design_file, header=T, sep = "\t", na.strings = "0", check.names=F)
rawcount = read.csv(rawcount_file, header=T, sep ="\t", check.names=F)

print(design)

name_sample= as.character(as.vector(design[,1]))
countMatrix = rawcount[,3:ncol(rawcount)]

# Iterate over each design

for (i in 2:ncol(design)) {
	
	name_folder = paste(out_path,names(design[i]),sep="/")
       
	 # Create output directory  
        
	if (!file.exists(name_folder)) { 
		system(paste("mkdir",name_folder,sep=" "))
	}

        current_design=design[,i]
        subsampleN=name_sample[!(is.na(current_design))]
        group = as.character(current_design)[!(is.na(current_design))]
        groupN = unique(group)
        current_countMatrix = NULL
        for (j in 1:length(subsampleN)) {
                current_countMatrix=cbind(current_countMatrix,countMatrix[,is.element(colnames(countMatrix),subsampleN[j])])
        }
        colnames(current_countMatrix)=subsampleN
        rownames(current_countMatrix)=rawcount[,1]
        libSize <- colSums(current_countMatrix)
        #geneSymbol=rawcount[,2]
	
        cat("Processing for the design\n")
        cat(paste("Name folder: ",name_folder,"\n",sep=""))
        cat(paste("Design : ",paste(subsampleN, group,sep="=",collapse=" ; "),"\n",spe=""))

	# Perform gene differential expressoin

       	perform_dge(current_countMatrix, group, count_limit, name_folder)
}



 
