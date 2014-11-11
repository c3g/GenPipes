# TODO:
# Bisulfite stuff...
# Small RNA
# Variation
# mappability



module load mugqic/R/3.0.2
module load mugqic/bwa/0.7.9a 
module load mugqic/samtools/0.1.19-gpfs
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.108
module load mugqic/tools/1.9
module load mugqic/bowtie2/2.2.2



##
######### hg19

## Fetch from iGenomes
GENOME="oryCun2"
SPECIE="Oryctolagus_cuniculus"
ROOT="$MUGQIC_INSTALL_HOME/genomes/${SPECIE}/${GENOME}"
mkdir -p $ROOT/fasta/bowtie2 $ROOT/fasta/bwa $ROOT/fasta/byChro $ROOT/annotation/mappability $ROOT/annotation/genes
cd $ROOT/fasta
wget ftp://ftp.ensembl.org/pub/release-75/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.75.dna_sm.toplevel.fa.gz
gunzip ${GENOME}.fa.gz 
ln -s ${GENOME}.fa ${GENOME}.fa.fa

##build chromosome files, dictionnaries, indexes, etc...
WG2ChromosomeFasta.awk ${GENOME}.fa byChro
samtools faidx ${GENOME}.fa
ln -s ${GENOME}.fa.fai ${GENOME}.fa.fa.fai
java -jar ${PICARD_HOME}/CreateSequenceDictionary.jar REFERENCE=${GENOME}.fa OUTPUT=${GENOME}.dict GENOME_ASSEMBLY=${GENOME}
ln -s ${GENOME}.dict ${GENOME}.dict.dict
ln -s ${GENOME}.dict ${GENOME}.fa.dict
ln -s ${GENOME}.dict ${GENOME}.fa.fa.dict
## Chrom size file, simple tab-delimied 1st col chr , 2nd col chr length. (already in iGenomes Annotation/Genes/ChromInfo.txt)
cut -f1,2 ${GENOME}.fa.fai >  ${GENOME}.chromsize.txt

# bowtie-build Index reference (already done for iGenomes)
cd $ROOT/fasta/bowtie2
ln -s ../${GENOME}.fa
ln -s ../${GENOME}.fa.fa
ln -s ../${GENOME}.fa.fai
ln -s ../${GENOME}.fa.fa.fai
ln -s ../${GENOME}.dict
bowtie2-build ${GENOME}.fa ${GENOME}.fa


# bwa index  reference (already done for iGenomes)
cd $ROOT/fasta/bwa
ln -s ../${GENOME}.fa
ln -s ../${GENOME}.fa.fai
ln -s ../${GENOME}.dict
bwa index -a bwtsw ${GENOME}.fa 


###Get annotations
#genes from ensembl
cd $ROOT/annotation/genes
wget ftp://ftp.ensembl.org/pub/release-75/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.75.gtf.gz
gunzip Oryctolagus_cuniculus.OryCun2.0.75.gtf.gz

# Gene size file for the gene model + gene id to symbol mapping. Easier to use GenomicRanges 

R --vanilla <<'EOF'
	#require(org.Hs.eg.db) # ORGANISM SPECIFIC!!!!
	#org.prefix = 'org.Hs.eg'	
	require(stringr)
	require(rtracklayer)
	options(stringsAsFactors=FALSE)
	x = import.gff('Oryctolagus_cuniculus.OryCun2.0.75.gtf')
	x$geneid = gsub('((gene_id |;)|\")','', str_extract(as.character(x$group),'gene_id.*?;') ) # MIGHT DIFFER FOR NON-UCSC GENOMES

	# Write the gene_id to symbol file
	g2s  = data.frame('ensembl'=unique(x$geneid),'geneSymbol'=unique(x$geneid))
	write.table(g2s,file='genes_geneid2Symbol.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

	# Prepare the gene lengths file
	gl = x
	gl = gl[gl$type=='exon',]
	gl =  as(gl,'GRanges') 		
	gl = split(gl,f=factor(gl$geneid))	 # split for GRanges signature turns into GRangesList!
	gl = reduce(gl) # reduce eliminates overlaps
	gl = sum ( width(gl) )
	gl = as.data.frame(gl)
	write.table(gl,file='genes_lengths.txt',col.names=FALSE,row.names=TRUE,sep='\t',quote=FALSE)



EOF


#rRNA
cd $ROOT/annotation/
# wget ftp://ftp.ensembl.org/pub/release-75/fasta/oryctolagus_cuniculus/ncrna/Oryctolagus_cuniculus.OryCun2.0.75.ncrna.fa.gz
# gunzip Oryctolagus_cuniculus.OryCun2.0.75.ncrna.fa.gz
awk ' BEGIN {pr=0} {if (substr($0,0,1) == ">") {x=split($0,na,"rRNA"); if (x > 1) {pr=1} else {pr=0}}; if (pr == 1) {print $0}} ' Oryctolagus_cuniculus.OryCun2.0.75.ncrna.fa > OryCun2.0.75.rRNA.fa
bwa index OryCun2.0.75.rRNA.fa


###GO file manually created using the better buny website:
# http://cptweb.cpt.wayne.edu/BB/index.php
# Please cite: Craig DB, Kannan S, Dombkowski AA. Augmented annotation and orthologue analysis for Oryctolagus cuniculus: Better Bunny. BMC Bioinformatics. 2012 May 8;13:84. PubMed PMID: 22568790; PubMed Central PMCID: PMC3424829.
	# dowloaded as bb.csv file in  $ROOT/annotation/genes
	# awk ' BEGIN {FS="\""} NR > 1 {x=split($0,go,"GO:") ; if (x > 1) {for (i=2;i<=x;i++) {y=split(go[i],goNum,")");print $2"\tGO:"goNum[1]}}} ' bb.csv  > OryCun2_genes2go.tsv


### Create  rRNA, tRNA, mtRNA mask for cufflinks (paste me to hg19)

# Manual: UCSC, rpsmk with filter repClass tRNA, rRNA
# Manual get all chrM from UCSCgenes table
# cat
#cat hg19_rpmsk_rRNA_tRNA.gtf chrM.gtf > hg19_rRNA_tRNA_chrM.gtf




# OR Alternative is MAthieu's:
# module load mugqic/tools
# gtf2geneSize.awk <...>









