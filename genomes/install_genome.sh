#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

init_install() {
  # '$MUGQIC_INSTALL_HOME_DEV' for development, '$MUGQIC_INSTALL_HOME' for production
  INSTALL_HOME=$MUGQIC_INSTALL_HOME_DEV

  INSTALL_DIR=$INSTALL_HOME/genomes/species/$SPECIES.$ASSEMBLY
  SOURCE_DIR=$INSTALL_DIR/source
  LOG_DIR=$INSTALL_DIR/log
  TIMESTAMP=`date +%FT%H.%M.%S`

  # Standard names for our genome files whatever the source
  GENOME_DIR=$INSTALL_DIR/genome
  GENOME_FASTA=$SPECIES.$ASSEMBLY.fa
  ANNOTATIONS_DIR=$INSTALL_DIR/annotations
  GTF=$SPECIES.$ASSEMBLY.$RELEASE.gtf
  NCRNA=$SPECIES.$ASSEMBLY.$RELEASE.ncrna.fa
  VCF=$SPECIES.$ASSEMBLY.$RELEASE.vcf.gz

  echo Installing genome for:
  echo species: $SPECIES
  echo assembly: $ASSEMBLY
  echo source: $SOURCE
  echo in: $INSTALL_DIR
  echo

  # Create install directory with permissions if necessary
  if [[ ! -d $INSTALL_DIR ]] ; then mkdir -p $INSTALL_DIR ; fi

  # Create subdirectories
  mkdir -p $SOURCE_DIR $LOG_DIR $GENOME_DIR $ANNOTATIONS_DIR
}

download_dir() {
  URL=$1

  # Download directory is the URL domain name
  echo $SOURCE_DIR/`echo $URL | perl -pe 's/^\S+:\/\/([^\/]*).*$/\1/'`
}

download_url() {
  URL=$1

  echo
  echo Downloading $URL...
  echo

  DOWNLOAD_DIR=`download_dir $URL`
  mkdir -p $DOWNLOAD_DIR
  cd $DOWNLOAD_DIR
  wget $URL
}

set_urls() {
  #
  # Ensembl (vertebrate species)
  #
  if [[ $SOURCE == "Ensembl" ]]
  then
    URL_PREFIX=ftp://ftp.ensembl.org/pub/release-$RELEASE
    GENOME_URL=$URL_PREFIX/fasta/${SPECIES,,}/dna/$SPECIES.$ASSEMBLY.dna.primary_assembly.fa.gz
    NCRNA_URL=$URL_PREFIX/fasta/${SPECIES,,}/ncrna/$SPECIES.$ASSEMBLY.ncrna.fa.gz
    GTF_URL=$URL_PREFIX/gtf/${SPECIES,,}/$SPECIES.$ASSEMBLY.$RELEASE.gtf.gz
    VCF_URL=$URL_PREFIX/variation/vcf/${SPECIES,,}/$SPECIES.vcf.gz
  
    # Before Ensembl release 76, release number was added in genome and ncrna file names
    if [ $RELEASE -lt 76 ]
    then
      GENOME_URL=${GENOME_URL/$SPECIES.$ASSEMBLY/$SPECIES.$ASSEMBLY.$RELEASE}
      NCRNA_URL=${NCRNA_URL/$SPECIES.$ASSEMBLY/$SPECIES.$ASSEMBLY.$RELEASE}
    fi
  
    # Check if a genome primary assembly is available for this species, otherwise use the toplevel assembly
    set +e
    wget --spider $GENOME_URL
    EXIT_CODE=$?
    set -e
    if [ $EXIT_CODE != 0 ]
    then
      echo Primary assembly not available for $SPECIES, use toplevel assembly instead
      GENOME_URL=${GENOME_URL/primary_assembly/toplevel}
    fi

    # Check if a VCF is available for this species
    set +e
    wget --spider $VCF_URL
    EXIT_CODE=$?
    set -e
    if [ $EXIT_CODE != 0 ]
    then
      echo VCF not available for $SPECIES
      VCF_URL=
    fi
  
  #
  # Ensembl Genomes (non-vertebrate species)
  #
  elif [[ $SOURCE == "Ensembl_Genomes" ]]
  then
    RELEASE_URL=ftp://ftp.ensemblgenomes.org/pub/release-$RELEASE
  
    # Retrieve Ensembl Genomes species information
    download_url $RELEASE_URL/species.txt
    cd `download_dir $RELEASE_URL/species.txt`
    SPECIES_LINE="`awk -F"\t" -v species=${SPECIES,,} '$2 == species' species.txt`"

    # Retrieve species division (Bacteria|Fungi|Metazoa|Plants|Protists)
    DIVISION=`echo "$SPECIES_LINE" | cut -f3 | sed "s/^Ensembl//"`
  
    # Escherichia coli bacteria file paths are different
    # Retrieve which bacteria collection it belongs to and adjust paths accordingly
    CORE_DB_PREFIX=`echo "$SPECIES_LINE" | cut -f13 | perl -pe "s/_core_${RELEASE}_\d+_\d+//"`
    if [[ $CORE_DB_PREFIX != ${SPECIES,,} ]]
    then
      EG_SPECIES=$CORE_DB_PREFIX/${SPECIES,,}
      EG_BASENAME=$SPECIES.`echo "$SPECIES_LINE" | cut -f6`.$RELEASE
    else
      EG_SPECIES=${SPECIES,,}
      EG_BASENAME=$SPECIES.$ASSEMBLY.$RELEASE
    fi
  
    URL_PREFIX=$RELEASE_URL/${DIVISION,}
    GENOME_URL=$URL_PREFIX/fasta/$EG_SPECIES/dna/$EG_BASENAME.dna.genome.fa.gz
    NCRNA_URL=$URL_PREFIX/fasta/$EG_SPECIES/ncrna/$EG_BASENAME.ncrna.fa.gz
    GTF_URL=$URL_PREFIX/gtf/$EG_SPECIES/$EG_BASENAME.gtf.gz
    if [[ `echo "$SPECIES_LINE" | cut -f13` == "Y"  ]]
    then
      VCF_URL=$URL_PREFIX/vcf/${SPECIES,,}/${SPECIES,,}.vcf.gz
    else
      echo VCF not available for $SPECIES
    fi
  fi
}

download_urls() {
  download_url $GENOME_URL
  download_url $GTF_URL
  download_url $NCRNA_URL
  if [ ! -z "${VCF_URL:-}" ]
  then
    download_url $VCF_URL
  fi
}

cmd_or_job() {
  CMD=$1

  # If genome is too big, run command in a separate job since login node memory is limited
  if [ `stat --printf="%s" $GENOME_DIR/$GENOME_FASTA` -gt 200000000 ]
  then
    echo
    echo Submitting $CMD as job...
    echo
    echo "${!CMD}" | qsub -m ae -M $JOB_MAIL -W umask=0002 -d $INSTALL_DIR -j oe -o $LOG_DIR/${CMD}_$TIMESTAMP.log -N $CMD.$GENOME_FASTA -l walltime=24:00:0 -q sw -l nodes=1:ppn=2
  else
    echo
    echo Running $CMD...
    echo
    echo "${!CMD}" | bash
  fi
}

create_picard_index() {
  echo
  echo Creating genome Picard sequence dictionary...
  echo
  module load mugqic/picard/1.108 mugqic/java
  java -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=$GENOME_DIR/$GENOME_FASTA OUTPUT=$GENOME_DIR/${GENOME_FASTA/.fa/.dict} GENOME_ASSEMBLY=${GENOME_FASTA/.fa} > $LOG_DIR/picard_$TIMESTAMP.log 2>&1
}

create_samtools_index() {
  echo
  echo Creating genome SAMtools FASTA index...
  echo
  module load mugqic/samtools/0.1.19
  samtools faidx $GENOME_DIR/$GENOME_FASTA > $LOG_DIR/samtools_$TIMESTAMP.log 2>&1
}

create_bwa_index() {
  echo
  echo Creating genome BWA index...
  echo
  BWA_CMD="\
INDEX_DIR=$GENOME_DIR/bwa_index && \
mkdir -p \$INDEX_DIR && \
ln -s -f -t \$INDEX_DIR ../$GENOME_FASTA
module load mugqic/bwa/0.7.10 && \
bwa index \$INDEX_DIR/$GENOME_FASTA > $LOG_DIR/bwa_$TIMESTAMP.log 2>&1"
  cmd_or_job BWA_CMD
}

create_bowtie2_tophat_index() {
  echo
  echo Creating genome Bowtie 2 index and gtf TopHat index...
  echo
  BOWTIE2_TOPHAT_CMD="\
INDEX_DIR=$GENOME_DIR/bowtie2_index && \
mkdir -p \$INDEX_DIR && \
ln -s -f -t \$INDEX_DIR ../$GENOME_FASTA
module load mugqic/bowtie/2.1.0 && \
bowtie2-build \$INDEX_DIR/$GENOME_FASTA \$INDEX_DIR/${GENOME_FASTA/.fa} > $LOG_DIR/bowtie2_$TIMESTAMP.log 2>&1 && \
INDEX_DIR=$INSTALL_DIR/annotations/gtf_tophat_index && \
mkdir -p \$INDEX_DIR && \
ln -s -f -t \$INDEX_DIR ../$GTF
module load mugqic/samtools/0.1.19 mugqic/tophat/2.0.11 && \
tophat --output-dir \$INDEX_DIR/tophat_out --GTF \$INDEX_DIR/$GTF --transcriptome-index=\$INDEX_DIR/${GTF/.gtf} $GENOME_DIR/bowtie2_index/${GENOME_FASTA/.fa} > $LOG_DIR/gtf_tophat_$TIMESTAMP.log 2>&1"
  cmd_or_job BOWTIE2_TOPHAT_CMD
}

create_ncrna_bwa_index() {
  echo
  echo Creating ncRNA BWA index...
  echo
  INDEX_DIR=$INSTALL_DIR/annotations/ncrna_bwa_index
  mkdir -p $INDEX_DIR
  ln -s -f -t $INDEX_DIR ../$NCRNA
  module load mugqic/bwa/0.7.10
  bwa index $INDEX_DIR/$NCRNA > $LOG_DIR/ncrna_bwa_$TIMESTAMP.log 2>&1
}

create_gene_annotations() {
  # Create gene annotation from GTF
  cd $INSTALL_DIR/annotations
  module load mugqic_dev/R/3.1.1
  R --no-restore --no-save<<EOF
suppressPackageStartupMessages(library(gqSeqUtils))
gtf.fn     = "$GTF"
annotation = "$SPECIES.$ASSEMBLY.$RELEASE"

genes = extractGeneAnnotationFromGtf(gtf.fn, feature.types = c("exon", "CDS"))
# NOTE: feature.types passed down to gqSeqUtils::calculateGeneLengthsFromGtf

write.table(genes[,c("featureID","gene.length"),],
  file = paste0(annotation, ".genes.length.tsv"), sep='\t', quote=F, row.names=F, col.names=F) # NOTE: this should not be necessary, all is in genes anyway

write.table(data.frame("ensembl"=genes[["featureID"]], "geneSymbol"=genes[["SYMBOL"]]),
  file = paste0(annotation, ".geneid2Symbol.tsv"), sep='\t', quote=F, row.names=F, col.names=F) # NOTE: this should not be necessary, all is in genes anyway

write.table(genes, file = paste0(annotation, ".genes.tsv"), sep='\t', quote=F, row.names=F, col.names=T)

EOF
}

get_vcf_dbsnp() {
  # Try to retrieve VCF dbSNP version (set +e since zgrep exit code != 0 if not found)
  set +e
  DBSNP=`zgrep -m1 -Po "##INFO=<ID=dbSNP_\d+" $ANNOTATIONS_DIR/$VCF | cut -f3 -d=`
  set -e
  echo Found VCF dbSNP version: $DBSNP
  if [ $DBSNP ]
  then
    # Add dbSNP version to VCF filename
    mv $ANNOTATIONS_DIR/$VCF $ANNOTATIONS_DIR/${VCF/.vcf.gz/.$DBSNP.vcf.gz}
  fi
}

copy_files() {
  # Uncompress files
  cd `download_dir $GENOME_URL`
  gunzip -c `basename $GENOME_URL` > $GENOME_DIR/$GENOME_FASTA
  cd `download_dir $GTF_URL`
  gunzip -c `basename $GTF_URL` > $ANNOTATIONS_DIR/$GTF
  cd `download_dir $NCRNA_URL`
  gunzip -c `basename $NCRNA_URL` > $ANNOTATIONS_DIR/$NCRNA

  if [ ! -z "${VCF_URL:-}" ]
  then
    cd `download_dir $VCF_URL`
    cp `basename $VCF_URL` $ANNOTATIONS_DIR/$VCF
    get_vcf_dbsnp
  fi
}

build_files() {
  # Create indexes
  create_picard_index
  create_samtools_index
  create_bwa_index
  create_bowtie2_tophat_index
  create_ncrna_bwa_index

  create_gene_annotations
}

create_genome_ini_file() {
  echo "\
[info]
scientific_name=$SPECIES
common_name=$COMMON_NAME
assembly=$ASSEMBLY
assembly_synonyms=$ASSEMBLY_SYNONYMS" > $INSTALL_DIR/genome.ini
}

install_genome() {

  SPECIES="$1"
  COMMON_NAME="$2"
  ASSEMBLY="$3"
  ASSEMBLY_SYNONYMS="$4"
  SOURCE="$5"
  RELEASE="$6"

  init_install
  set_urls
  download_urls
  copy_files
  build_files
  create_genome_ini_file

  # Add permissions
  chmod ug+rwX,o+rX $INSTALL_DIR
}
