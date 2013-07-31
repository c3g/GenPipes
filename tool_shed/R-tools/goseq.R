# Perform a GO term analysis with package goseq
# Maxime Caron Nov 2011
# Borrowed some argument code from Mathieu Bourgey
# usage : Rscript goseq.R -d path_dge_results_file -c comma_delimited_columns -g organism -k known_reference -o output_dir
# modified  by Mathieu Bourgey Feb 2013
# modified  by Mathieu Bourgey July 2013


usage=function(errM) {
        cat("\nUsage : Rscript goseq.R [option] <Value>\n")
        cat("       -d        : file containing the dge results\n")
        cat("       -c        : columns to use for input (first column = geneSymbol list, second column= adjusted p-values\n")
        cat("       -t        : annotation reference name (default org.Hs.eg.db)\n")
        cat("       -s        : specie UCSC reference name (default hg19)\n")
	cat("       -k        : known genome reference file\n")
	cat("       -m        : maximum number of pathway return (default all)\n")
	cat("       -v        : DGE method 0: edger/deseq ; 1: cuffdiff (default 0)\n")
	cat("       -i        : gene Name ID (default geneSymbol)\n")
	cat("       -G        : optional: Nonnative Go file path (default NULL)\n")
	cat("       -a        : optional: gene length file path (default NULL)\n")
	cat("       -o        : output directory\n")
        cat("       -h        : this help\n\n")
        stop(errM)
}


##################################
ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 10) {
usage("missing arguments")
}

## defult arg values
file=""
columns=""
organism=""
out_path=""
annotation="org.Hs.eg.db"
maxP= -1
method=0
specie="hg19"
gene_path=NULL
go_path=""
gene_type="geneSymbol"
useNNgo=FALSE

## get arg variables
for (i in 1:length(ARG)) {
        if (ARG[i] == "-d") {
                file=ARG[i+1]
        } else if (ARG[i] == "-c") {
                columns=ARG[i+1]
        } else if (ARG[i] == "-t") {
                annotation=ARG[i+1]
	} else if (ARG[i] == "-s") {
                specie=ARG[i+1]
	} else if (ARG[i] == "-o") {
                out_path=ARG[i+1]
	} else if (ARG[i] == "-k") {
                known_ref=ARG[i+1]
	} else if (ARG[i] == "-m") {
                maxP=as.numeric(ARG[i+1])
	} else if (ARG[i] == "-v") {
                method=as.numeric(ARG[i+1])
	} else if (ARG[i] == "-G") {
                go_path=ARG[i+1]
	} else if (ARG[i] == "-a") {
                gene_path=ARG[i+1]
	} else if (ARG[i] == "-i") {
                gene_type=ARG[i+1]
	} else if (ARG[i] == "-h") {
                usage("")
        }
}
## check arg consitency
if (!(file.exists(file))) {
	stop("Input file not found") 
}
if (out_path == "") {
	stop("Output directory not found")
}


library('goseq')

testAnno=find.package(annotation,quiet=T)
print(testAnno)
repo = grep(".*\\..*\\.LENGTH", as.data.frame(data(package = "geneLenDataBase")$results, stringsAsFactors = FALSE)$Item, ignore.case = TRUE, value = TRUE)
repo = matrix(unlist(strsplit(repo, "\\.")), ncol = 3, byrow = TRUE)
valid_ids = sapply(split(repo[, 2], repo[, 1]), paste, collapse = ",")
supportedSpecies=names(valid_ids)
supportedGeneIDList=strsplit(valid_ids,",")

if (length(testAnno) && (specie %in% supportedSpecies)) {
	require(annotation,character.only=T)
	if (gene_type %in% unlist(supportedGeneIDList[supportedSpecies %in% specie])) {
		print(paste("annotation file:",annotation))
		print(paste("specie:",specie))
		print(paste("gene ID:",gene_type))
	} else if (file.exists(go_path) && file.exists(gene_path)) {
		print(paste("Using Non-native Gene Identifier",gene_path,"and category test",go_path))
		useNNgo=TRUE
	} else {
		stop("Wrong gene ID while no usable Non-native file are provided")
	}
} else if (file.exists(go_path) && file.exists(gene_path)) {
	print(paste("Using Non-native Gene Identifier",gene_path,"and category test",go_path))
	useNNgo=TRUE
} else {
stop("Wrong annotation file or specie, while no usable Non-native file are provided")
}

set.seed(123456789)

d1<-read.table(file, header=T, sep="\t", quote="", stringsAsFactors=F)
toUse<-unlist(strsplit(columns,","))
selecT<-c(as.numeric(toUse[1]),as.numeric(toUse[2]))
d2<-d1[,selecT]

if(method == 1) {
kgX<-read.table(known_ref, header=T, sep="\t", quote="", stringsAsFactors=F)
tmp<-merge(d2,kgX,by.x=1,by.y=1)
d2<-tmp[,c(3,2)]
}
#data_externe[order(data_externe[,5]),]
d2<-d2[order(d2[,2]),]

head(d2)
is.significant<-function(x) ifelse(x<0.05,1,0)
if(sum(is.significant(d2[,2])==1) == 0) {
stop("No significant adjusted p-values found")
}
d3<-cbind(d2[,1], is.significant(d2[,2]))
de<-subset(d3,d3[,2]==1)
gene.vector = as.integer(unique(d3[,1]) %in% unique(de[,1]))
names(gene.vector) = unique(d3[,1])
if (useNNgo) {
	goTable=read.table(go_path,header=F)
	geneLenPre=read.table(gene_path,header=F)
	geneTable=unfactor(geneLenPre[,2])
	names(geneTable)=unfactor(geneLenPre[,1])
	pwf = nullp(gene.vector, bias.data=geneTable)
	GO.wall =  goseq(pwf,gene2cat = goTable)
} else {
	pwf = nullp(gene.vector, specie, gene_type)
	GO.wall = goseq(pwf, specie, gene_type)
}
head(GO.wall)
enriched.GO = cbind(GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < 0.05], GO.wall$over_represented_pvalue[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < 0.05])
head(enriched.GO)
library(GO.db)
if(dim(enriched.GO)[1] == 0) {
	print("No significant FDR < 0.05 enriched categories")
	q("no",0)
} else {
if (maxP == -1) maxP=dim(enriched.GO)[1]
}
#write.table("Gene Ontology Analysis", out_path, row.names=F, col.names=F, quote=F)
write.table(paste("Enriched category","FDR < 0.05 filtered p-value","GOID","Term","Ontology","Definition","Synonym", sep="\t"), out_path, append=F, row.names=F, col.names=F, quote=F)
for (i in 1:maxP) {
f<-GOTERM[[enriched.GO[i,1]]]
write.table(paste(i, enriched.GO[i,2], GOID(f), Term(f), Ontology(f), Definition(f), Synonym(f)[i], sep="\t"), out_path, append=T, row.names=F, col.names=F, quote=F)
}

#getgo(d3[1:1000,15], "mm9", "geneSymbol", fetch.cats=c("GO:CC","GO:BP","GO:MF"))
