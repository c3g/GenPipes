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
SPECIES=Saccharomyces_cerevisiae
ASSEMBLY=Source_assembly_version0.0  ## TO BE MODIFIED WITH e.g. GRCh37, CanFam3.1, TAIR10, etc.
ASSEMBLY=R64-1-1

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV  ## TO BE MODIFIED IF NECESSARY

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/genomes/species/$SPECIES/$ASSEMBLY

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir -p $INSTALL_DIR
  chmod ug+rwX,o+rX $INSTALL_DIR
fi

#SOURCE=(Ensembl|UCSC)  ## TO BE MODIFIED WITH SPECIFIC SOURCE
SOURCE=Ensembl

#
# Ensembl
#
if [[ $SOURCE == "Ensembl" ]]
then

ENSEMBL_RELEASE=75  ## TO BE MODIFIED IF NECESSARY
ENSEMBL_URL_PREFIX=ftp://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE

DNA_URL=$ENSEMBL_URL_PREFIX/fasta/${SPECIES,,}/dna/$SPECIES.$ASSEMBLY.$ENSEMBL_RELEASE.dna.toplevel.fa.gz
GTF_URL=$ENSEMBL_URL_PREFIX/gtf/${SPECIES,,}/$SPECIES.$ASSEMBLY.$ENSEMBL_RELEASE.gtf.gz
rRNA_URL=$ENSEMBL_URL_PREFIX/fasta/${SPECIES,,}/ncrna/$SPECIES.$ASSEMBLY.$ENSEMBL_RELEASE.ncrna.fa.gz

fi

cd $INSTALL_DIR

echo Installing genome assembly...
mkdir -p $INSTALL_DIR/dna
cd $INSTALL_DIR/dna
wget $DNA_URL
gunzip `basename $DNA_URL`
ln -s `basename ${DNA_URL%.gz}` $SPECIES.$ASSEMBLY.fa

echo Creating genome SAM dictionary...
module load mugqic/picard mugqic/java
java -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=$SPECIES.$ASSEMBLY.fa OUTPUT=$SPECIES.$ASSEMBLY.dict GENOME_ASSEMBLY=$SPECIES.$ASSEMBLY

echo Creating genome FASTA index...
module load mugqic/samtools
samtools faidx $SPECIES.$ASSEMBLY.fa

echo Creating genome bwa index...
mkdir -p $INSTALL_DIR/dna/bwa
cd $INSTALL_DIR/dna/bwa
ln -s ../$SPECIES.$ASSEMBLY.fa
module load mugqic/bwa
bwa index $SPECIES.$ASSEMBLY.fa
