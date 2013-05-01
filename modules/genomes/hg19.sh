# TODO:
# Bisulfite stuff...
# Small RNA
# Variation
# mappability



module load mugqic/R
module load mugqic/bwa
module load mugqic/samtools
module load mugqic/picard




##
######### hg19

## Fetch from iGenomes
GENOME="hg19"
cd $MUGQIC_INSTALL_HOME/genomes
wget "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz"
tar -xvf Homo_sapiens_UCSC_hg19.tar.gz
ROOT="$MUGQIC_INSTALL_HOME/genomes/Homo_sapiens/UCSC/hg19"

# bowtie-build Index reference (already done for iGenomes)
#    bowtie-build /lb/project/mugqic/epigenome/genome_files/hg19/fasta/hg19.fasta /lb/project/mugqic/epigenome/genome_files/hg19/fasta/hg19.fasta # example Maxime
# bwa index  reference (already done for iGenomes)
#   bwa index -a bwtsw fasta/hg19.fasta # Example Maxime

# faidx reference (already done for iGenomes)
#   samtools faidx fasta/hg19.fasta # example Maxime

# reference dictionary (already done for iGenomes)
#cd $ROOT/Sequence/WholeGenomeFasta
#java -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=genome.fa  


## Link WholeGenomeFasta/* in aligner indices  dirs
for dir in $ROOT/Sequence/{BWA,Bowtie,Bowtie2}Index
do
 cd $dir
 ln -s "../WholeGenomeFasta/genome.fa.fai"
 ln -s "../WholeGenomeFasta/genome.dict"
done




## Chrom size file, simple tab-delimied 1st col chr , 2nd col chr length. (already in iGenomes Annotation/Genes/ChromInfo.txt)
# OR cut -f1,2 $ROOT/Sequence/WholeGenomeFasta/genome.fa.fai > $ROOT/Sequence/WholeGenomeFasta/genome.chromsize.txt
# OR  sed 1d fasta/hg19.fasta.dict | awk '{print $2"\t"$3}' | sed -e 's/SN://g' | sed -e 's/LN://g' > hg19.chromsize.txt # example Maxime

## rRNA sequence + BWA index them (either from the Broad OR from iGenomes) (ln -s to). Put into $ROOT/Sequence/AbundantSequences/rRNA.fa
# cat $ROOT/Sequence/AbundantSequences/hum5SrDNA.fa $ROOT/Sequence/AbundantSequences/humRibosomal.fa > $ROOT/Sequence/AbundantSequences/rRNA.fa
# OR :
wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/examples/rRNA/rRNA.tar.gz
tar -zxvf rRNA.tar.gz 
rm human_all_rRNA.fasta.*
mv human_all_rRNA.fasta $ROOT/Sequence/AbundantSequences/rRNA.fa
bwa index $ROOT/Sequence/AbundantSequences/rRNA.fa

## Gene models aka gtf.   (Already done for iGenomes, Annotation/Genes/genes.gtf)
# OR:
# too complicated... ha

## BioC TxDB
# not working... use the txdb from Bioc???


## Gene size file for the gene model + gene id to symbol mapping. Easier to use GenomicRanges 
cd $ROOT/Annotation/Genes
R --vanilla <<'EOF'
	require(stringr)
	require(rtracklayer)
	options(stringsAsFactors=FALSE)
	x = import.gff('genes.gtf')
	x$geneid = gsub('((gene_id |;)|\")','', str_extract(as.character(x$group),'gene_id.*?;') )
	g2s  = data.frame('ensembl'=unique(x$geneid),'geneSymbol'=unique(x$geneid))
	write.table(g2s,file='genes_geneid2Symbol.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)
	x = x[x$type=='exon',]
	x =  as(x,'GRanges') 		
	x = split(x,f=factor(x$geneid))	
	x = reduce(x) # reduce eliminates overlaps
	x = sum ( width(x) )
	x = as.data.frame(x)
	write.table(x,file='genes_lengths.txt',col.names=FALSE,row.names=TRUE,sep='\t',quote=FALSE)
EOF
cd -
# OR Alternative is MAthieu's:
# module load mugqic/tools
# gtf2geneSize.awk <...>









