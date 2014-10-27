
##
#### Genome Installation on abacus (copy me to install script when finished)

# source: http://ciliate.org/index.php/home/downloads
# http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5911&lvl=3&lin=f&keep=1&srchmode=1&unlock

## Setup 
ROOT=/sb/programs/analyste/genomes/Tetrahymena_thermophila/


## Download source files
mkdir -p $ROOT/source && cd source
echo "http://www.ciliate.org/system/downloads" > "README.txt"
wget http://www.ciliate.org/system/downloads/T_thermophila_dec2011.gff3
wget http://www.ciliate.org/system/downloads/T_thermophila_oct2008_proteins.fasta
wget http://www.ciliate.org/system/downloads/T_thermophila_oct2008_CDS.fasta
wget http://www.ciliate.org/system/downloads/T_thermophila_oct2008_gene.fasta
wget http://www.ciliate.org/system/downloads/T_thermophila_oct2008_assembly.fasta
wget http://www.ciliate.org/system/downloads/T_thermophila_domain_info.txt
wget http://www.ciliate.org/system/downloads/T_thermophila_final_gene_attributes.txt

## .dict, .fai file and all that
module load mugqic
cd $ROOT
mkdir -p WholeGenomeFasta && cd WholeGenomeFasta
ln -s ../source/T_thermophila_oct2008_assembly.fasta genome.fa
# Create reference dictionary

module load mugqic/picard/1.96;java -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=genome.fa OUTPUT=genome.dict
sed 1d WholeGenomeFasta/genome.dict | awk '{print $2"\t"$3}' | sed -e 's/SN://g' | sed -e 's/LN://g' > genome.chromsize.txt
# .fai

module load mugqic/samtools/0.1.18;samtools faidx genome.fa

## BWA

cd $ROOT
mkdir -p bwa && cd bwa
ln -s ../source/T_thermophila_oct2008_assembly.fasta genome.fa
module load mugqic/bwa/0.7.4;bwa index genome.fa
module load mugqic/samtools/0.1.18;samtools faidx genome.fa


## Bowtie2 index

module load mugqic/bowtie/2.1.0
cd $ROOT
mkdir -p Bowtie2Index && cd Bowtie2Index
ln -s ../source/T_thermophila_oct2008_assembly.fasta genome.fa
bowtie2-build genome.fa genome


## gtf

cd $ROOT
mkdir -p annotation && cd annotation
ln -s ../source/T_thermophila_dec2011.gff3
module load mugqic/cufflinks;gffread -E T_thermophila_dec2011.gff3 -T -o- > T_thermophila_dec2011.gtf

## create referenceEnsemble2symbol

awk ' BEGIN {FS="\""} {print $4 "\t" $6} ' T_thermophila_dec2011.gtf | sort -u > T_thermophila_dec2011_geneID2geneName.tsv

## Ribosomal fasta

cd $ROOT

mkdir -p annotations/Ribosomal_fasta && cd annotations/Ribosomal_fasta
ln -s ../../source/T_thermophila_dec2011.gff3

module load mugqic/hmmer/2.3.2;module load mugqic/rnammer/1.2;rnammer -S euk -m lsu,ssu,tsu -xml seq.xml -gff seq.gff -h seq.hmmreport -f ribosomal_Tthermophila.fasta < T_thermophila_oct2008_assembly.fasta
module load mugqic/bwa/0.7.4;bwa index ribosomal_Tthermophila.fasta



## Go annotations

cd $ROOT

mkdir -p annotations/GO_annotation && cd annotations/GO_annotation
ln -s ../../source/tthermophila.obo.txt
## This file was send by e-mail, the window caractere need to be take off.
awk ' {x=length($5) ; goID="GO:" ; for (i=1;i<=(7-x);i++) {goID=goID "0"} ; print $2 "\t" goID $5 } ' tthermophila.obo.txt > tthermophila2GO.tsv













