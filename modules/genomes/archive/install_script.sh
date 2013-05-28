## The purpose of this script is to provide a template for standardized installation of new genomes.
## iGenomes already contains many useful files such as fastas, bwa/bowtie indices, fai, dict, etc.
## We could try to follow the iGenome standards for file names , e.g. gtf as Annotation/Genes/gene.gtf
##
##
##
##
module load mugqic/R
module load mugqic/bwa
module load mugqic/samtools
module load mugqic/picard




##
######### hg19

## Fetch from iGenomes
GENOME="hg19"
cd $MUGQIC_INSTALL_HOME/genomes
#wget "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz"
#tar -xvf Homo_sapiens_UCSC_hg19.tar.gz
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

## Link aligner indices in WholeGenome
cd $ROOT/Sequence/WholeGenomeFasta
# BWA
for fn in `ls  ../BWAIndex/genome.fa.*`
do
ln -sf $fn
done
# Bowtie2
for fn in `ls  ../Bowtie2Index/genome.*.*`
do
ln -sf $fn
done
# Bowtie
for fn in `ls  ../BowtieIndex/genome.*.*`
do
ln -sf $fn
done
cd -

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




# Bisulfite stuff...
# Small RNA
# Variation
# mappability











##
######### canFam2
GENOME="canFam2"

## Fetch from iGenomes
cd $MUGQIC_INSTALL_HOME/genomes
#wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Canis_familiaris/UCSC/canFam2/Canis_familiaris_UCSC_canFam2.tar.gz
#tar -xvf Canis_familiaris_UCSC_canFam2.tar.gz
ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/UCSC/canFam2"


## Link aligner indices in WholeGenome
cd $ROOT/Sequence/WholeGenomeFasta
# BWA
for fn in `ls  ../BWAIndex/genome.fa.*`
do
ln -sf $fn
done
# Bowtie2
for fn in `ls  ../Bowtie2Index/genome.*.*`
do
ln -sf $fn
done
# Bowtie
for fn in `ls  ../BowtieIndex/genome.*.*`
do
ln -sf $fn
done
cd -


## rRNA sequence + BWA index them (either from the Broad OR from iGenomes) (ln -s to). Put into $ROOT/Sequence/AbundantSequences/rRNA.fa
# http://en.wikipedia.org/wiki/RRNA#Eukaryotes
# TODO: BiomaRt???




##
######### canFam3
GENOME="canFam3"

## Fetch from iGenomes
cd $MUGQIC_INSTALL_HOME/genomes
ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Sequence/Archives"
mkdir -p $ROOT
cd $ROOT
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.1.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.2.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.3.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.4.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.5.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.6.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.7.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.8.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.9.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.10.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.11.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.12.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.13.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.14.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.15.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.16.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.17.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.18.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.19.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.20.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.21.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.22.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.23.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.24.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.25.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.26.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.27.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.28.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.29.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.31.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.32.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.33.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.34.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.35.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.36.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.37.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.38.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.X.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.chromosome.MT.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.71.dna.nonchromosomal.fa.gz

ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Sequence/Chromosomes/"
mkdir -p $ROOT

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 
do
  zcat Canis_familiaris.CanFam3.1.71.dna.chromosome.$i.fa.gz > $ROOT/Chr$i.fa
done
zcat Canis_familiaris.CanFam3.1.71.dna.nonchromosomal.fa.gz > $ROOT/Chromosomes/ChrUn.fa
zcat Canis_familiaris.CanFam3.1.71.dna.chromosomal.MT.fa.gz > $ROOT/Chromosomes/ChrM.fa
zcat Canis_familiaris.CanFam3.1.71.dna.chromosomal.X.fa.gz > $ROOT/Chromosomes/ChrX.fa
cd $ROOT

ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Sequence/WholeGenomeFasta/"
mkdir -p $ROOT
cat * >> $ROOT/genome.fa
sed 's/\.1/_1/g' $ROOT/genome.fa > $ROOT/genome.new.fa
mv $ROOT/genome.fa $ROOT/genome.old.fa
mv $ROOT/genome.new.fa $ROOT/genome.fa
module add mugqic/samtools/0.1.19
samtools faidx $ROOT/genome.fa

cd $ROOT
module add mugqic/java/oracle-jdk1.7.0_15
module add mugqic/picard/1.88
java -Xmx8g -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=$ROOT/genome.fa O=$ROOT/genome.dict

cd $P_INDEX
P_INDEX=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Sequence/BWAIndex/"
mkdir -p $P_INDEX
ln -s $ROOT/genome.fa $P_INDEX/genome.fa
module add mugqic/bwa/0.7.4
bwa index $P_INDEX/genome.fa
ln -s $P_INDEX/genome.fa.* $ROOT

cd $P_INDEX
P_INDEX=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Sequence/Bowtie2Index/"
mkdir -p $P_INDEX
ln -s $ROOT/genome.fa $P_INDEX/genome.fa
module add mugqic/bowtie/2.0.6
bowtie2-build $P_INDEX/genome.fa $P_INDEX/genome
ln -s $P_INDEX/genome.fa.* $ROOT

ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Annotation/rRNA/"
mkdir -p $ROOT
cp $MUGQIC_INSTALL_HOME/genomes/Canis_familiaris/UCSC/canFam2Annotation/rRNA/* Annotation/rRNA

ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Annotation/Genes/"
mkdir -p $ROOT
cd $ROOT
wget ftp://ftp.ensembl.org/pub/release-71/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.71.gtf.gz
zcat Canis_familiaris.CanFam3.1.71.gtf.gz > genes.gtf
rm Canis_familiaris.CanFam3.1.71.gtf.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/database/chromInfo.txt.gz
zcat chromInfo.txt.gz > ChromInfo.txt
rm chromInfo.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/database/refFlat.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/database/refGene.txt.gz
zcat refGene.txt.gz > refGene.txt
rm refGene.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/database/refSeqSummary.txt.gz
zcat refSeqSummary.txt.gz > refSeqSummary.txt
rm refSeqSummary.txt.gz

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


ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Annotation/Variation/gvf/"
mkdir -p $ROOT
cd $ROOT
wget ftp://ftp.ensembl.org/pub/release-71/variation/gvf/canis_familiaris/Canis_familiaris.gvf.gz
wget ftp://ftp.ensembl.org/pub/release-71/variation/gvf/canis_familiaris/Canis_familiaris_failed.gvf.gz
wget ftp://ftp.ensembl.org/pub/release-71/variation/gvf/canis_familiaris/Canis_familiaris_incl_consequences.gvf.gz
wget ftp://ftp.ensembl.org/pub/release-71/variation/gvf/canis_familiaris/Canis_familiaris_structural_variations.gvf.gz
wget ftp://ftp.ensembl.org/pub/release-71/variation/gvf/canis_familiaris/README
ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Annotation/Variation/vcf/"
mkdir -p $ROOT
cd $ROOT
wget ftp://ftp.ensembl.org/pub/release-71/variation/vcf/canis_familiaris/Canis_familiaris.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-71/variation/vcf/canis_familiaris/Canis_familiaris_incl_consequences.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-71/variation/vcf/canis_familiaris/Canis_familiaris_structural_variations.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-71/variation/vcf/canis_familiaris/README


ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Annotation/cDNA/"
mkdir -p $ROOT
cd $ROOT
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/cdna/Canis_familiaris.CanFam3.1.71.cdna.abinitio.fa.gz
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/cdna/Canis_familiaris.CanFam3.1.71.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/cdna/README

ROOT=$MUGQIC_INSTALL_HOME"/genomes/Canis_familiaris/canFam3/Annotation/ncRNA/"
mkdir -p $ROOT
cd $ROOT
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/cdna/Canis_familiaris.CanFam3.1.71.ncrna.fa.gz
wget ftp://ftp.ensembl.org/pub/release-71/fasta/canis_familiaris/cdna/README




## rRNA sequence + BWA index them (either from the Broad OR from iGenomes) (ln -s to). Put into $ROOT/Sequence/AbundantSequences/rRNA.fa
# http://en.wikipedia.org/wiki/RRNA#Eukaryotes
# TODO: BiomaRt???


R --vanilla <<'EOF'
	require(stringr)
	require(rtracklayer)
	options(stringsAsFactors=FALSE)
	x = import.gff('Schizosaccharomyces_pombe.ASM294v1.18.gtf')
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
