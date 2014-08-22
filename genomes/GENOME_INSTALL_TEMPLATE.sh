#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a genome install script template which should be copied and used for
# consistency between genome paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


SPECIES=Species_scientific_name  ## TO BE MODIFIED WITH e.g. Homo_sapiens, Canis_familiaris, Arabidopsis_thaliana, etc.
#SPECIES=Saccharomyces_cerevisiae
#SPECIES=Caenorhabditis_elegans
#SPECIES=Mus_musculus
SPECIES=Homo_sapiens
#SPECIES=Escherichia_coli_str_k_12_substr_dh10b
#SPECIES=Arabidopsis_thaliana
#SPECIES=Schizosaccharomyces_pombe
ASSEMBLY=Source_assembly_version0.0  ## TO BE MODIFIED WITH e.g. GRCh37, CanFam3.1, TAIR10, etc.
#ASSEMBLY=R64-1-1
#ASSEMBLY=WBcel235
#ASSEMBLY=GRCm38
ASSEMBLY=GRCh38
#ASSEMBLY=ASM1942v1
#ASSEMBLY=TAIR10
#ASSEMBLY=ASM294v2
#SOURCE=(Ensembl|Ensembl_Genomes)  ## TO BE MODIFIED WITH SPECIFIC SOURCE
#SOURCE=Ensembl_Genomes
SOURCE=Ensembl
RELEASE=76  ## TO BE MODIFIED WITH SPECIFIC Ensembl(_Genomes) RELEASE IF NECESSARY
BASENAME=$SPECIES.$ASSEMBLY.$RELEASE

# '$MUGQIC_INSTALL_HOME_DEV' for development, '$MUGQIC_INSTALL_HOME' for production
INSTALL_HOME=$MUGQIC_INSTALL_HOME_DEV  ## TO BE MODIFIED IF NECESSARY

INSTALL_DIR=$INSTALL_HOME/genomes/species/$SPECIES.$ASSEMBLY
SOURCE_DIR=$INSTALL_DIR/source/$SOURCE
LOG_DIR=$INSTALL_DIR/log
TIMESTAMP=`date +%FT%H.%M.%S`

# Standard name for our genome FASTA file whatever the source
GENOME_DIR=$INSTALL_DIR/genome
GENOME_BASENAME=$SPECIES.$ASSEMBLY
GENOME_FASTA=$GENOME_BASENAME.fa
GTF=$BASENAME.gtf
NCRNA=$BASENAME.ncrna.fa

echo Installing genome for:
echo species: $SPECIES
echo assembly: $ASSEMBLY
echo source: $SOURCE
echo in: $INSTALL_DIR
echo

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]] ; then mkdir -p $INSTALL_DIR ; fi

# Create subdirectories
mkdir -p $SOURCE_DIR $LOG_DIR $GENOME_DIR $INSTALL_DIR/annotations

cd $SOURCE_DIR

#
# Ensembl (vertebrate species)
#
if [[ $SOURCE == "Ensembl" ]]
then

URL_PREFIX=ftp://ftp.ensembl.org/pub/release-$RELEASE
#GENOME_URL=$URL_PREFIX/fasta/${SPECIES,,}/dna/$GENOME_BASENAME.dna.(primary_assembly|toplevel).fa.gz  ## TO BE MODIFIED WITH SPECIFIC GENOME URL
GENOME_URL=$URL_PREFIX/fasta/${SPECIES,,}/dna/$GENOME_BASENAME.dna.primary_assembly.fa.gz
NCRNA_URL=$URL_PREFIX/fasta/${SPECIES,,}/ncrna/$GENOME_BASENAME.ncrna.fa.gz
GTF_URL=$URL_PREFIX/gtf/${SPECIES,,}/$BASENAME.gtf.gz

#
# Ensembl Genomes (non-vertebrate species)
#
elif [[ $SOURCE == "Ensembl_Genomes" ]]
then

RELEASE_URL=ftp://ftp.ensemblgenomes.org/pub/release-$RELEASE
SPECIES_URL=$RELEASE_URL/species.txt
wget $SPECIES_URL
# Retrieve species division (Bacteria|Fungi|Metazoa|Plants|Protists)
DIVISION=`cut -f2,3 species.txt | grep -i ${SPECIES,,} | cut -f2 | sed "s/^Ensembl//"`

# Escherichia coli bacteria file paths are different
# Retrieve which bacteria collection it belongs to
CORE_DB_PREFIX=`cut -f2,13 species.txt | grep -i ${SPECIES,,} | cut -f2 | perl -pe "s/_core_${RELEASE}_\d+_\d+//"`
if [[ $CORE_DB_PREFIX != ${SPECIES,,} ]]
then
EG_SPECIES=$CORE_DB_PREFIX/${SPECIES,,}
EG_BASENAME=$SPECIES.`cut -f2,6 species.txt | grep -i ${SPECIES,,} | cut -f2`.$RELEASE
else
EG_SPECIES=${SPECIES,,}
EG_BASENAME=$BASENAME
fi

URL_PREFIX=$RELEASE_URL/${DIVISION,}
GENOME_URL=$URL_PREFIX/fasta/$EG_SPECIES/dna/$EG_BASENAME.dna.genome.fa.gz
NCRNA_URL=$URL_PREFIX/fasta/$EG_SPECIES/ncrna/$EG_BASENAME.ncrna.fa.gz
GTF_URL=$URL_PREFIX/gtf/$EG_SPECIES/$EG_BASENAME.gtf.gz
fi

echo
echo Downloading $GENOME_URL...
echo
#wget $GENOME_URL
#gunzip -c `basename $GENOME_URL` > $GENOME_DIR/$GENOME_FASTA

echo
echo Downloading $GTF_URL...
echo
#wget $GTF_URL
#gunzip -c `basename $GTF_URL` > $INSTALL_DIR/annotations/$GTF

echo
echo Downloading $NCRNA_URL...
echo
wget $NCRNA_URL
gunzip -c `basename $NCRNA_URL` > $INSTALL_DIR/annotations/$NCRNA

# Create indexes directories
for INDEX_DIR in bowtie2_index bwa_index picard_index sam_index
do
  mkdir -p $GENOME_DIR/$INDEX_DIR
  cd $GENOME_DIR/$INDEX_DIR
  ln -s -f -t $GENOME_DIR/$INDEX_DIR ../$GENOME_FASTA
done

echo
echo Creating genome Picard sequence dictionary...
echo
INDEX_DIR=$GENOME_DIR/picard_index
module load mugqic/picard/1.108 mugqic/java
java -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=$INDEX_DIR/$GENOME_FASTA OUTPUT=$INDEX_DIR/$GENOME_BASENAME.dict GENOME_ASSEMBLY=$GENOME_BASENAME > $LOG_DIR/picard_$TIMESTAMP.log 2>&1

echo
echo Creating genome SAMtools FASTA index...
echo
INDEX_DIR=$GENOME_DIR/sam_index
module load mugqic/samtools/0.1.19
samtools faidx $INDEX_DIR/$GENOME_FASTA > $LOG_DIR/samtools_$TIMESTAMP.log 2>&1

echo
echo Creating ncRNA BWA index...
echo
INDEX_DIR=$INSTALL_DIR/annotations/ncrna_bwa_index
mkdir -p $INDEX_DIR
ln -s -f -t $INDEX_DIR ../$NCRNA
module load mugqic/bwa/0.7.10
bwa index $INDEX_DIR/$NCRNA > $LOG_DIR/ncrna_bwa_$TIMESTAMP.log 2>&1

BWA_CMD="INDEX_DIR=$GENOME_DIR/bwa_index && \
module load mugqic/bwa/0.7.10 && \
bwa index \$INDEX_DIR/$GENOME_FASTA > $LOG_DIR/bwa_$TIMESTAMP.log 2>&1"

BOWTIE2_TOPHAT_CMD="
INDEX_DIR=$GENOME_DIR/bowtie2_index && \
module load mugqic/bowtie/2.1.0 && \
bowtie2-build \$INDEX_DIR/$GENOME_FASTA \$INDEX_DIR/$GENOME_BASENAME > $LOG_DIR/bowtie2_$TIMESTAMP.log 2>&1 && \
INDEX_DIR=$INSTALL_DIR/annotations/gtf_tophat_index && \
mkdir -p \$INDEX_DIR && \
ln -s -f -t \$INDEX_DIR ../$GTF
module load mugqic/samtools/0.1.19 mugqic/tophat/2.0.11 && \
tophat --output-dir \$INDEX_DIR/tophat_out --GTF \$INDEX_DIR/$GTF --transcriptome-index=\$INDEX_DIR/$BASENAME $GENOME_DIR/bowtie2_index/$GENOME_BASENAME > $LOG_DIR/gtf_tophat_$TIMESTAMP.log 2>&1"

# If genome is too big, create index in a separate job since login node memory is limited
if [ `stat --printf="%s" $GENOME_DIR/$GENOME_FASTA` -gt 200000000 ]
then
QSUB="qsub -m ae -M $JOB_MAIL -W umask=0002 -d $INSTALL_DIR -j oe -o $LOG_DIR/__JOB_$TIMESTAMP.log -N __JOB.$GENOME_BASENAME -l walltime=24:00:0 -q sw -l nodes=1:ppn=2"
fi

for CMD in BWA_CMD BOWTIE2_TOPHAT_CMD
do
# If QSUB variable is not defined, run command directly
if [ -z "${QSUB:-}" ]
then
  echo
  echo Running ${CMD}...
  echo
  echo "${!CMD}" | bash
else
  echo
  echo Submitting ${CMD} as job...
  echo
  echo "${!CMD}" | ${QSUB//__JOB/$CMD}
fi
done

# Add permissions
chmod ug+rwX,o+rX $INSTALL_DIR
