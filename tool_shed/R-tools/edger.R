# Performs differential gene expression with EdgeR
# Written by Maxime Caron, Mathieu Bourgey & Greg Voisin - Dec 2011
# Usage : Rscript edger.R -d path_design -c path_rawcountfile -o output_dir

library(edgeR)

# Usage

usage=function(errM) {
	cat("\nUsage : Rscript edger.R [option] <Value>\n")
	cat("       -d        : design file\n")
	cat("       -c        : raw count file\n")
	cat("       -o        : output directory\n")
	cat("       -h        : this help\n\n")
	stop(errM)
}


set.seed(123456789)
perform_dge=function(d, count_limit=1, path, numb, genSymbol) {

# Retain row which have > count_limit

d<-d[rowSums(d$counts) > count_limit,]

# Normalize read counts

dTMM<-calcNormFactors(d)
de.common=NULL
#print(dTMM[1:5,])

if(numb == 1) {
	dTMM2=estimateCommonDisp(dTMM)
	dTMM3<-estimateTagwiseDisp(dTMM2)
	de.common<-exactTest(dTMM3)
}
else if(numb == 2) {
	dTMM<-estimateGLMCommonDisp(dTMM,matrix(1,2,1),method="deviance",robust=TRUE)
	de.common<-exactTest(dTMM)
}

print(de.common)

summary_TMM=summary(decideTestsDGE(de.common))
top_TMM=topTags(de.common, n=dim(d)[1], sort.by = "logFC")
if (!is.null(de.common$genes)) {
	rownames(top_TMM$table)=rownames(de.common$table)[match(top_TMM$table$genes,de.common$genes[,1])]
}
if(numb == 1) {
count_order_pseudo=dTMM3$pseudo.counts[match(rownames(top_TMM$table), rownames(d$count)),]
colnames(count_order_pseudo)=paste(colnames(count_order_pseudo),"norm",sep=".")
}
count_order_real=d$count[match(rownames(top_TMM$table), rownames(d$count)),]
colnames(count_order_real)=paste(colnames(count_order_real),"raw",sep=".")
id=rownames(top_TMM$table)
geneN=genSymbol[names(genSymbol) %in% id]
colnames(top_TMM$table)=c("gene_symbol", "log_conc","log_fc", "edger.p-value","edger.adj.p-value")
top_TMM$table[,2:3]=round(top_TMM$table[,2:3], digits=3)
top_TMM$table[,4]=as.numeric(format(top_TMM$table[,4], digits=2))
top_TMM$table[,5]=as.numeric(format(top_TMM$table[,5], digits=2))
top_TMM$table[,1]=geneN[match(id,names(geneN))]
if(numb == 1) {
data_externe=cbind(id,top_TMM$table, round(count_order_pseudo, digits=0), round(count_order_real, digits=0))
}
else {
data_externe=cbind(id,top_TMM$table, round(count_order_real, digits=0))
}
write.table(data_externe[order(data_externe[,5]),], paste(path,"edger_results.csv",sep="/") , quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)


}


##################################

ARG = commandArgs(trailingOnly = T)
## default arg values
count_limit=9
fpath="."
design_file=""
rawcount_file=""
genlist_file=""
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
# ## check arg consitency
# if (!(file.exists(design_file))) {
# 	usage("Error : Design file not found")
# }
# if (!(file.exists(rawcount_file))) {
# 	usage("Error : Raw count file not found")
# }
# if (!(file.exists(genlist_file))) {
# 	usage("Error : gene file not found")
# }
# if (out_path == "") {
# 	usage("Error : Output directory not specified")
# }
# tmpFP=strsplit(fpath,"")
# if (tmpFP[[1]][length(tmpFP[[1]])] == "/" ) {
# 	bckS=""
# } else {
# 	bckS="/"
# }
# tmpOP=strsplit(out_path,"")
# if (tmpOP[[1]][length(tmpOP[[1]])] == "/" ) {
# 	out_path=paste(tmpOP[[1]][1:(length(tmpOP[[1]]-1))],collapse="")
# }

design = read.csv2(design_file, header=T, sep = "\t", na.strings = "0", check.names=F)
rawcount = read.csv(rawcount_file, header=T, sep ="\t", check.names=F)
genL=cbind(rawcount[,1:2])

print(design)

name_sample= as.character(as.vector(design[,1]))
countMatrix = rawcount[,3:ncol(rawcount)]

# Iterate over each design

for (i in 2:ncol(design)) {
	
	name_folder = paste(out_path,names(design[i]),sep="/")
       
	 # Create output directory  
        
	if (!file.exists(name_folder)) { 
		dir.create(name_folder, showWarnings=F, recursive=T)
	}

        current_design=design[,i]
        subsampleN=name_sample[!(is.na(current_design))]
        group = as.character(current_design)[!(is.na(current_design))]
        groupN = unique(group)
        current_countMatrix_tmp = NULL
        for (j in 1:length(subsampleN)) {
                current_countMatrix_tmp=cbind(current_countMatrix_tmp,countMatrix[,is.element(colnames(countMatrix),subsampleN[j])])
        }
	current_countMatrix=current_countMatrix_tmp[rowSums(current_countMatrix_tmp) > 0,]
        colnames(current_countMatrix)=subsampleN
        rownames(current_countMatrix)=rawcount[rowSums(current_countMatrix_tmp) > 0,1]
        libSize <- colSums(current_countMatrix)
        geneSymbol=genL[genL[,1] %in% rownames(current_countMatrix),2]
	geneID=genL[genL[,1] %in% rownames(current_countMatrix),1]
	names(geneSymbol)=geneID
	
        cat("Processing for the design\n")
        cat(paste("Name folder: ",name_folder,"\n",sep=""))
        cat(paste("Design : ",paste(subsampleN, group,sep="=",collapse=" ; "),"\n",spe=""))

	# Both groups have replicates

	dge.list=DGEList(counts=current_countMatrix, group=group , lib.size=libSize, genes=geneID)

	if((table(group)[1] > 1) & (table(group)[2] > 1)) {
        	perform_dge(dge.list, count_limit=length(subsampleN), name_folder, 1, geneSymbol)
        } 
	
	# Both groups have no replicates

	else if((table(group)[1] == 1) & (table(group)[2] == 1)) {
        	perform_dge(dge.list, count_limit=length(subsampleN),name_folder, 2, geneSymbol)
	}

	# Group 1 has replicates but not group 2	

	else if((table(group)[1] > 1) & (table(group)[2] == 1)) {
                perform_dge(dge.list, count_limit=length(subsampleN),name_folder, 1, geneSymbol)
	}

	# Group 2 has replicates but not group 1
	
	else if((table(group)[1] == 1) & (table(group)[2] > 1)) {
                perform_dge(dge.list, count_limit=length(subsampleN),name_folder, 1, geneSymbol)
        }
}
