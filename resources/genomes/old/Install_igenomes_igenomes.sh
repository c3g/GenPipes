#Install a genome from i-genomes

#Example for Rattus norvegicus

## Set your variables

##

genomeName="Rattus_norvegicus";
build="rn5";
bwaVersion="6";
bowtieVersion="2";
bowtieDir=`echo "bowtie"$bowtieVersion`
bwaDir=`echo "bwa"$bwaVersion`
INSTALL_PATH=$MUGQIC_INSTALL_HOME/genomes/$genomeName
dictName=`echo $genomeName".dict"`

##

cd $INSTALL_PATH

## change the path for your genomes and the gunzip name and tar name

wget --ftp-user=igenome --ftp-password=G3nom3s4u ftp://ftp.illumina.com/Rattus_norvegicus/UCSC/rn5/Rattus_norvegicus_UCSC_rn5.tar.gz
gunzip Rattus_norvegicus_UCSC_rn5.tar.gz
tar -xvf Rattus_norvegicus_UCSC_rn5.tar

mv $genomeName/UCSC .
mv UCSC source
mkdir -p fasta
mkdir -p annotations
cd $INSTALL_PATH/fasta
cp ../source/$build/Sequence/WholeGenomeFasta/genome.fa .
mv genome.fa $genomeName.fa


mkdir -p $bowtieDir
mkdir -p $bwaDir

module load mugqic/tools mugqic/bwa mugqic/bowtie mugqic/samtools mugqic/picard

# Index reference with samtools faidx reference (already done for iGenomes)
##samtools faidx $genomeName.fa
cp ../source/rn5/Sequence/WholeGenomeFasta/genome.fa.fai .
mv genome.fa.fai $genomeName".fai"

# reference dictionary (already done for iGenomes)
module load mugqic/java mugqic/picard && java -jar ${PICARD_HOME}/CreateSequenceDictionary.jar REFERENCE=$genomeName.fa OUTPUT=$dictName


# Generate tab delimited chromosome size file (already done for iGenomes)

cp ../source/$build/GenomeStudio/$genomeName/UCSC-"$build"/ChromInfo.txt .
mv ChromInfo.txt chromosome_size.txt

cd $INSTALL_PATH/fasta/$bowtieDir

cp ../../source/$build/Sequence/Bowtie2Index/genome.* .

cd $INSTALL_PATH/fasta/$bwaDir

cp ../../source/$build/Sequence/BWAIndex/genome.* .

#Annotations

cd $INSTALL_PATH/annotations
cp ../source/$build/Annotation/Genes/genes.gtf .
module load mugqic/tools/1.6 && gtf2geneSize.awk genes.gtf genes.tsv
module load mugqic/tools/1.6 && gtf2tmpMatrix.awk genes.gtf gene2geneName.tsv











