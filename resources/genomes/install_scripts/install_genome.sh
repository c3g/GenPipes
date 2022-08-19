#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

module_bismark=mugqic/bismark/0.21.0
module_bowtie=mugqic/bowtie/1.2.2
module_bowtie2=mugqic/bowtie2/2.3.5
module_bwa=mugqic/bwa/0.7.17
module_java=mugqic/java/openjdk-jdk1.8.0_72
module_mugqic_tools=mugqic/mugqic_tools/2.4.0
module_mugqic_R_packages=mugqic/mugqic_R_packages/1.0.6
module_perl=mugqic/perl/5.22.1
module_picard=mugqic/picard/2.0.1
module_python=mugqic/python/2.7.14
module_R=mugqic/R_Bioconductor/3.6.1_3.10
module_samtools=mugqic/samtools/1.11
module_star=mugqic/star/2.7.8a
module_tabix=mugqic/tabix/0.2.6
module_tophat=mugqic/tophat/2.0.14
module_ucsc=mugqic/ucsc/v387
module_hicup=mugqic/HiCUP/v0.7.2
module_kallisto=mugqic/kallisto/0.44.0

HOST=`hostname`

# Ensure to use 'grep' from CVMFS to avoid errors caused by different grep versions
grep_cvmfs=/cvmfs/soft.mugqic/yum/centos7/1.0/usr/bin/grep

init_install() {
  # '$MUGQIC_INSTALL_HOME_DEV' for development, '$MUGQIC_INSTALL_HOME' for production
  if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
fi

  INSTALL_DIR=${!INSTALL_HOME}/genomes/species/$SPECIES.$ASSEMBLY
  DOWNLOAD_DIR=$INSTALL_DIR/downloads
  LOG_DIR=$INSTALL_DIR/log
  TIMESTAMP=`date +%FT%H.%M.%S`

  # Standard names for our genome files whatever the source
  GENOME_DIR=$INSTALL_DIR/genome
  GENOME_FASTA=$SPECIES.$ASSEMBLY.fa
  ANNOTATIONS_DIR=$INSTALL_DIR/annotations
  GTF=$SPECIES.$ASSEMBLY.$SOURCE$VERSION.gtf
  NCRNA=$SPECIES.$ASSEMBLY.$SOURCE$VERSION.ncrna.fa
  RRNA=$SPECIES.$ASSEMBLY.$SOURCE$VERSION.rrna.fa
  VCF=$SPECIES.$ASSEMBLY.$SOURCE$VERSION.vcf.gz
  GO=$SPECIES.$ASSEMBLY.$SOURCE$VERSION.GO.tsv
  CDNA=$SPECIES.$ASSEMBLY.$SOURCE$VERSION.cdna.fa

  echo "Installing genome for:"
  echo "species: $SPECIES"
  echo "assembly: $ASSEMBLY"
  echo "source: $SOURCE"
  echo "in: $INSTALL_DIR"
  echo

  # Create install directory with permissions if necessary
  if [[ ! -d $INSTALL_DIR ]] ; then mkdir -p $INSTALL_DIR ; fi

  # Create subdirectories
  mkdir -p $DOWNLOAD_DIR $LOG_DIR $GENOME_DIR $ANNOTATIONS_DIR
}

is_url_valid() {
  URL=$1

  set +e
  wget `dirname $URL`/ -O - | grep `basename $URL`
  EXIT_CODE=$?
  set -e

  return $EXIT_CODE
}

download_path() {
  URL=$1

  # Download path is the URL without protocol prefix
  echo $DOWNLOAD_DIR/`echo $URL | perl -pe 's,^[^:]*://,,'`
}

download_url() {
  URL=$1

  echo
  echo "Downloading $URL..."
  echo

  cd $DOWNLOAD_DIR
  wget -r -c $URL
}

set_urls() {
  #
  # Ensembl (vertebrate species)
  #
  if [[ $SOURCE == "Ensembl" ]]
  then
    URL_PREFIX=ftp://ftp.ensembl.org/pub/release-$VERSION
    GENOME_URL=$URL_PREFIX/fasta/${SPECIES,,}/dna/$SPECIES.$ASSEMBLY.dna.primary_assembly.fa.gz
    NCRNA_URL=$URL_PREFIX/fasta/${SPECIES,,}/ncrna/$SPECIES.$ASSEMBLY.ncrna.fa.gz
    CDNA_URL=$URL_PREFIX/fasta/${SPECIES,,}/cdna/$SPECIES.$ASSEMBLY.cdna.all.fa.gz
    GTF_URL=$URL_PREFIX/gtf/${SPECIES,,}/$SPECIES.$ASSEMBLY.$VERSION.gtf.gz
    VCF_URL=$URL_PREFIX/variation/vcf/${SPECIES,,}/$SPECIES.vcf.gz
    VCF_TBI_URL=$VCF_URL.tbi
    BIOMART_MART=ENSEMBL_MART_ENSEMBL
    # Retrieve species short name e.g. "mmusculus" for "Mus_musculus"
    if [[ $(echo "${SPECIES}" | awk -F"_" '{print NF-1}') == 2 ]]
    then
        SPECIES_2ND_NAME=`echo ${SPECIES#*_}`
        SPECIES_SHORT_NAME=`echo ${SPECIES:0:1}${SPECIES_2ND_NAME:0:1}${SPECIES##*_} | tr [:upper:] [:lower:]`
    else
        SPECIES_SHORT_NAME=`echo ${SPECIES:0:1}${SPECIES#*_} | tr [:upper:] [:lower:]`
    fi
    BIOMART_DATASET=${SPECIES_SHORT_NAME}_gene_ensembl
    BIOMART_GENE_ID=ensembl_gene_id
    BIOMART_GO_ID=go_id
    BIOMART_GO_NAME=name_1006
    BIOMART_GO_DEFINITION=definition_1006

    # Before Ensembl release 76, release number was added in genome and ncrna file names
    if [ $VERSION -lt 76 ]
    then
      GENOME_URL=${GENOME_URL/$SPECIES.$ASSEMBLY/$SPECIES.$ASSEMBLY.$VERSION}
      NCRNA_URL=${NCRNA_URL/$SPECIES.$ASSEMBLY/$SPECIES.$ASSEMBLY.$VERSION}
      CDNA_URL=${CDNA_URL/$SPECIES.$ASSEMBLY/$SPECIES.$ASSEMBLY.$VERSION}
    fi

    # Check if a genome primary assembly is available for this species, otherwise use the toplevel assembly
    if ! is_url_valid $GENOME_URL
    then
      echo "Primary assembly not available for $SPECIES, use toplevel assembly instead"
      GENOME_URL=${GENOME_URL/primary_assembly/toplevel}
    fi

    # Check if a VCF is available for this species
    if ! is_url_valid $VCF_URL
    then
      echo "VCF not available for $SPECIES"
      VCF_URL=
    fi

    # Check if a VCF tabix index is available for this species
    if ! is_url_valid $VCF_TBI_URL
    then
      echo "VCF tabix index not available for $SPECIES"
      VCF_TBI_URL=
    fi

  #
  # Ensembl Genomes (non-vertebrate species)
  #
  elif [[ $SOURCE == "EnsemblGenomes" ]]
  then
    RELEASE_URL=ftp://ftp.ensemblgenomes.org/pub/release-$VERSION

    # Retrieve Ensembl Genomes species information
    SPECIES_URL=$RELEASE_URL/species.txt
    download_url $SPECIES_URL
    SPECIES_LINE="`awk -F"\t" -v species=${SPECIES,,} '$2 == species' $(download_path $SPECIES_URL)`"

    # Retrieve species division (Bacteria|Fungi|Metazoa|Plants|Protists)
    DIVISION=`echo "$SPECIES_LINE" | cut -f3 | sed "s/^Ensembl//"`

    # Escherichia coli bacteria file paths are different
    # Retrieve which bacteria collection it belongs to and adjust paths accordingly
    CORE_DB_PREFIX=`echo "$SPECIES_LINE" | cut -f13 | perl -pe "s/_core_${VERSION}_\d+_\d+//"`
    if [[ $CORE_DB_PREFIX != ${SPECIES,,} ]]
    then
      EG_SPECIES=$CORE_DB_PREFIX/${SPECIES,,}
      EG_BASENAME=$SPECIES.`echo "$SPECIES_LINE" | cut -f6`.$VERSION
    else
      EG_SPECIES=${SPECIES,,}
      EG_BASENAME=$SPECIES.$ASSEMBLY.$VERSION
    fi

    URL_PREFIX=$RELEASE_URL/${DIVISION,}
    GENOME_URL=$URL_PREFIX/fasta/$EG_SPECIES/dna/$EG_BASENAME.dna.genome.fa.gz
    NCRNA_URL=$URL_PREFIX/fasta/$EG_SPECIES/ncrna/$EG_BASENAME.ncrna.fa.gz
    CDNA_URL=$URL_PREFIX/fasta/$EG_SPECIES/cdna/$EG_BASENAME.cdna.all.fa.gz
    GTF_URL=$URL_PREFIX/gtf/$EG_SPECIES/$EG_BASENAME.gtf.gz
    if [[ `echo "$SPECIES_LINE" | cut -f8` == "Y"  ]]
    then
      VCF_URL=$URL_PREFIX/vcf/${SPECIES,,}/${SPECIES,,}.vcf.gz
    else
      echo "VCF not available for $SPECIES"
    fi
    # Retrieve species short name e.g. "athaliana" for "Arabidopsis_thaliana"
    SPECIES_SHORT_NAME=`echo ${SPECIES:0:1}${SPECIES#*_} | tr [:upper:] [:lower:]`
    BIOMART_DATASET=${SPECIES_SHORT_NAME}_eg_gene
    BIOMART_GENE_ID=ensembl_gene_id
    BIOMART_GO_ID=go_id
    BIOMART_GO_NAME=name_1006
    BIOMART_GO_DEFINITION=definition_1006

  #
  # Ensembl Plants
  #
  elif [[ $SOURCE == "EnsemblPlants" ]]
  then
    RELEASE_URL=ftp://ftp.ensemblgenomes.org/pub/plants/release-$VERSION

    # Retrieve Ensembl Genomes species information
    SPECIES_URL=$RELEASE_URL/species_EnsemblPlants.txt
    download_url $SPECIES_URL
    SPECIES_LINE="`awk -F"\t" -v species=${SPECIES,,} '$2 == species' $(download_path $SPECIES_URL)`"


    CORE_DB_PREFIX=`echo "$SPECIES_LINE" | cut -f13 | perl -pe "s/_core_${VERSION}_\d+_\d+//"`
    if [[ $CORE_DB_PREFIX != ${SPECIES,,} ]]
    then
      EG_SPECIES=$CORE_DB_PREFIX/${SPECIES,,}
      EG_BASENAME=$SPECIES.`echo "$SPECIES_LINE" | cut -f6`.$VERSION
    else
      EG_SPECIES=${SPECIES,,}
      EG_BASENAME=$SPECIES.$ASSEMBLY
    fi

    URL_PREFIX=$RELEASE_URL
    GENOME_URL=$URL_PREFIX/fasta/$EG_SPECIES/dna/$EG_BASENAME.dna.toplevel.fa.gz
    NCRNA_URL=$URL_PREFIX/fasta/$EG_SPECIES/ncrna/$EG_BASENAME.ncrna.fa.gz
    CDNA_URL=$URL_PREFIX/fasta/$EG_SPECIES/cdna/$EG_BASENAME.cdna.all.fa.gz
    GTF_URL=$URL_PREFIX/gtf/$EG_SPECIES/$EG_BASENAME.$VERSION.gtf.gz
    if [[ `echo "$SPECIES_LINE" | cut -f8` == "Y"  ]]
    then
      VCF_URL=$URL_PREFIX/vcf/${SPECIES,,}/${SPECIES,,}.vcf.gz
    else
      echo "VCF not available for $SPECIES"
    fi
    # Retrieve species short name e.g. "athaliana" for "Arabidopsis_thaliana"
    SPECIES_SHORT_NAME=`echo ${SPECIES:0:1}${SPECIES#*_} | tr [:upper:] [:lower:]`
    BIOMART_DATASET=${SPECIES_SHORT_NAME}_eg_gene
    BIOMART_GENE_ID=ensembl_gene_id
    BIOMART_GO_ID=go_accession
    BIOMART_GO_NAME=go_name_1006
    BIOMART_GO_DEFINITION=go_definition_1006

  #
  # UCSC
  #
  elif [[ $SOURCE == "UCSC" ]]
  then
    GENOME_URL=http://hgdownload.soe.ucsc.edu/goldenPath/$ASSEMBLY/bigZips/$ASSEMBLY.2bit
  fi
}

download_urls() {
  download_url $GENOME_URL

  # Annotations are not downloaded for UCSC genomes
  if [[ $SOURCE != "UCSC" ]]
  then
    download_url $GTF_URL
    download_url $NCRNA_URL
    download_url $CDNA_URL
    if [ ! -z "${VCF_URL:-}" ]
    then
      download_url $VCF_URL
    fi
    if [ ! -z "${VCF_TBI_URL:-}" ]
    then
      download_url $VCF_TBI_URL
    fi
  fi
}

# Test genome size to decide whether indexing requires more cores or memory
is_genome_big() {
  GENOME_SIZE_LIMIT=2000000

  if [ `stat --printf="%s" $GENOME_DIR/$GENOME_FASTA` -gt $GENOME_SIZE_LIMIT ]
  then
    return 0
  else
    return 1
  fi
}

# Test if a list of files given as parameters exist, are regular files and are not zero size
is_up2date() {
  # By default, files are up to date : 0 is true, 1 is false
  IS_UP2DATE=0

  for f in $@
  do
    if [[ (! -f $f && ! -d $f) || ! -s $f ]]
    then
      IS_UP2DATE=1
    fi
  done

  return $IS_UP2DATE
}

cmd_or_job() {
  CMD=$1
  CORES=${2:-1}  # Nb cores = 2nd param if defined else 1
  JOB_PREFIX=${3:-$CMD}  # Job prefix = 3rd param if defined else cmd name
  JOB_NAME=$JOB_PREFIX
  MEMORY=${4:-31G}

  # If genome is too big, run command in a separate job since login node memory is limited
  if is_genome_big
  then
    echo
    echo "Submitting $JOB_PREFIX as job..."
    echo
    if [[ $HOST == abacus* ]]; then
      echo "${!CMD}" | qsub -m ae -M $JOB_MAIL -A $RAP_ID -W umask=0002 -d $INSTALL_DIR -j oe -o $LOG_DIR/${JOB_PREFIX}_$TIMESTAMP.log -N $JOB_NAME -q qfat256 -l pmem=256000m -l walltime=12:00:0 -l nodes=1:ppn=$CORES
    else      
      echo "#! /bin/bash 
      ${!CMD}" | \
      sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $INSTALL_DIR -o $LOG_DIR/${JOB_PREFIX}_$TIMESTAMP.log -J $JOB_NAME --time=12:00:0 --mem=$MEMORY -N 1 -c $CORES
    fi
  else
    echo
    echo "Running $JOB_PREFIX..."
    echo
    echo "${!CMD}" | bash
  fi
}
	
create_picard_index() {
  GENOME_DICT=$GENOME_DIR/${GENOME_FASTA/.fa/.dict}

  if ! is_up2date $GENOME_DICT
  then
    echo
    echo "Creating genome Picard sequence dictionary..."
    echo
    PICARD_CMD="\
module load $module_picard $module_java && \
java -jar \$PICARD_HOME/picard.jar CreateSequenceDictionary REFERENCE=$GENOME_DIR/$GENOME_FASTA OUTPUT=$GENOME_DICT GENOME_ASSEMBLY=${GENOME_FASTA/.fa} > $LOG_DIR/picard_$TIMESTAMP.log 2>&1"
    cmd_or_job PICARD_CMD 8
  else
    echo
    echo "Genome Picard sequence dictionary up to date... skipping"
    echo
  fi
}

create_samtools_index() {
  if ! is_up2date $GENOME_DIR/$GENOME_FASTA.fai
  then
    echo
    echo "Creating genome SAMtools FASTA index..."
    echo
    module load $module_samtools
    samtools faidx $GENOME_DIR/$GENOME_FASTA > $LOG_DIR/samtools_$TIMESTAMP.log 2>&1
    awk '{print $1\t$2}' $GENOME_DIR/$GENOME_FASTA.fai > $GENOME_DIR/$GENOME_FASTA.tsv
  else
    echo
    echo "Genome SAMtools FASTA index up to date... skipping"
    echo
  fi
}

create_bismark_genome_reference() {
  BISMARK_INDEX_DIR=$GENOME_DIR/bismark_index
  BISMARK_GENOME_DICT=$BISMARK_INDEX_DIR/${GENOME_FASTA/.fa/.dict}
  if ! is_up2date $BISMARK_INDEX_DIR/$GENOME_FASTA
  then
    echo
    echo "Creating custom Genome Reference for Bismark :"
    echo "Appending sequences of lambda phage and pUC19 vector to genome reference FASTA file..."
    echo
    mkdir -p $BISMARK_INDEX_DIR
    cat $GENOME_DIR/$GENOME_FASTA ${!INSTALL_HOME}/genomes/lambda_phage.fa ${!INSTALL_HOME}/genomes/pUC19.fa > $BISMARK_INDEX_DIR/$GENOME_FASTA
  else
    echo
    echo "Custom Genome Reference for Bismark is up to date... skipping"
    echo
  fi
  if ! is_up2date $BISMARK_INDEX_DIR/$GENOME_FASTA.fai
  then
    echo
    echo "Creating SAMtools FASTA index for the custom Genome Reference for Bismark..."
    echo "loading $module_samtools"
    module load $module_samtools
    echo "launching samtools faidx $BISMARK_INDEX_DIR/$GENOME_FASTA > $LOG_DIR/samtools_for_bismark_$TIMESTAMP.log 2>&1"
    samtools faidx $BISMARK_INDEX_DIR/$GENOME_FASTA > $LOG_DIR/samtools_for_bismark_$TIMESTAMP.log 2>&1 
    echo "done"
  else
    echo
    echo "Genome SAMtools FASTA index for the custom Genome Reference for Bismark is up to date... skipping"
    echo
  fi
  if ! is_up2date $BISMARK_INDEX_DIR/${GENOME_FASTA/.fa/.dict}
  then
    echo
    echo "Creating genome Picard sequence dictionary for the custom Genome Reference for Bismark..."
    echo
    INDEX_CMD="\
module load $module_picard $module_java && \
java -jar \$PICARD_HOME/picard.jar CreateSequenceDictionary REFERENCE=$BISMARK_INDEX_DIR/$GENOME_FASTA OUTPUT=$BISMARK_INDEX_DIR/${GENOME_FASTA/.fa/.dict} GENOME_ASSEMBLY=${GENOME_FASTA/.fa} > $LOG_DIR/picard_for_bismark_$TIMESTAMP.log 2>&1"
    cmd_or_job INDEX_CMD 8
  else
    echo
    echo "Genome Picard sequence dictionary for the custom Genome Reference for Bismark is up to date... skipping"
    echo
  fi
  if ! is_up2date $BISMARK_INDEX_DIR/Bisulfite_Genome
  then
    echo
    echo "Creating Bisulfite Genome Reference with Bismark..."
    echo
    BISMARK_CMD="\
module load $module_bismark $module_bowtie2 && \
LOG=$LOG_DIR/bismark_genome_preparation_$TIMESTAMP.log && \
ERR=$LOG_DIR/bismark_genome_preparation_$TIMESTAMP.err &&\
bismark_genome_preparation $BISMARK_INDEX_DIR > \$LOG 2> \$ERR && \
chmod -R ug+rwX,o+rX $BISMARK_INDEX_DIR \$LOG"
    cmd_or_job BISMARK_CMD 8
  else
    echo
    echo "Bismark Bisulfite Genome Reference is up to date... skipping"
    echo
  fi
  if ! is_up2date $ANNOTATIONS_DIR/${ASSEMBLY}_wgbs_bin100bp_GC.bed
  then
    echo
    echo "Creating GC content file for the Bisulfite Genome Reference"
    echo
  BINGC_CMD="\
module load $module_python $module_mugqic_tools && \
LOG=$LOG_DIR/wgbs_bin100bp_GC_$TIMESTAMP.log && \
\$PYTHON_TOOLS/getFastaBinedGC.py -s 100 -r $BISMARK_INDEX_DIR/$GENOME_FASTA -o $ANNOTATIONS_DIR/${ASSEMBLY}_wgbs_bin100bp_GC.bed > \$LOG 2>&1" # && \
#chmod -R ug+rwX,o+rX $BISMARK_INDEX_DIR \$LOG"
    cmd_or_job BINGC_CMD 8
  else
    echo
    echo "GC content file for Bisulfite Genome Reference is up to date... skipping"
    echo
  fi
}

create_bwa_index() {
  INDEX_DIR=$GENOME_DIR/bwa_index
  if ! is_up2date $INDEX_DIR/$GENOME_FASTA.sa
  then
    echo
    echo "Creating genome BWA index..."
    echo
    BWA_CMD="\
mkdir -p $INDEX_DIR && \
ln -s -f -t $INDEX_DIR ../$GENOME_FASTA && \
module load $module_bwa && \
LOG=$LOG_DIR/bwa_$TIMESTAMP.log && \
bwa index $INDEX_DIR/$GENOME_FASTA > \$LOG 2>&1 && \
chmod -R ug+rwX,o+rX $INDEX_DIR \$LOG"
    cmd_or_job BWA_CMD 6
  else
    echo
    echo "Genome BWA index up to date... skipping"
    echo
  fi
}

create_bowtie_index() {
  BOWTIE_INDEX_DIR=$GENOME_DIR/bowtie_index
  BOWTIE_INDEX_PREFIX=$BOWTIE_INDEX_DIR/${GENOME_FASTA/.fa}

  if ! is_up2date $BOWTIE_INDEX_PREFIX.[1-4].ebwt $BOWTIE_INDEX_PREFIX.rev.[12].ebwt
  then
    echo
    echo "Creating genome Bowtie index and gtf TopHat index..."
    echo
    BOWTIE_CMD="\
mkdir -p $BOWTIE_INDEX_DIR && \
ln -s -f -t $BOWTIE_INDEX_DIR ../$GENOME_FASTA && \
module load $module_bowtie && \
LOG=$LOG_DIR/bowtie_$TIMESTAMP.log && \
ERR=$LOG_DIR/bowtie_$TIMESTAMP.err && \
bowtie-build $BOWTIE_INDEX_DIR/$GENOME_FASTA $BOWTIE_INDEX_PREFIX > \$LOG 2> \$ERR" # && \
#chmod -R ug+rwX,o+rX $BOWTIE_INDEX_DIR \$LOG \$ERR"
  cmd_or_job BOWTIE_CMD 4
  else
    echo
    echo "Genome Bowtie index up to date... skipping"
    echo
  fi
}

create_bowtie2_tophat_index() {
  BOWTIE2_INDEX_DIR=$GENOME_DIR/bowtie2_index
  BOWTIE2_INDEX_PREFIX=$BOWTIE2_INDEX_DIR/${GENOME_FASTA/.fa}
  TOPHAT_INDEX_DIR=$ANNOTATIONS_DIR/gtf_tophat_index
  TOPHAT_INDEX_PREFIX=$TOPHAT_INDEX_DIR/${GTF/.gtf}

  if is_up2date $ANNOTATIONS_DIR/$GTF
  then
    if ! is_up2date $BOWTIE2_INDEX_PREFIX.[1-4].bt2 $BOWTIE2_INDEX_PREFIX.rev.[12].bt2 $TOPHAT_INDEX_PREFIX.[1-4].bt2 $TOPHAT_INDEX_PREFIX.rev.[12].bt2
    then
      echo
      echo "Creating genome Bowtie 2 index and gtf TopHat index..."
      echo
      BOWTIE2_TOPHAT_CMD="\
mkdir -p $BOWTIE2_INDEX_DIR && \
ln -s -f -t $BOWTIE2_INDEX_DIR ../$GENOME_FASTA && \
module load $module_bowtie2 && \
LOG=$LOG_DIR/bowtie2_$TIMESTAMP.log && \
ERR=$LOG_DIR/bowtie2_$TIMESTAMP.err && \
bowtie2-build $BOWTIE2_INDEX_DIR/$GENOME_FASTA $BOWTIE2_INDEX_PREFIX > \$LOG 2> \$ERR && \
chmod -R ug+rwX,o+rX $BOWTIE2_INDEX_DIR \$LOG \$ERR && \
mkdir -p $TOPHAT_INDEX_DIR && \
ln -s -f -t $TOPHAT_INDEX_DIR ../$GTF && \
module load $module_samtools $module_tophat && \
LOG=$LOG_DIR/gtf_tophat_$TIMESTAMP.log && \
ERR=$LOG_DIR/gtf_tophat_$TIMESTAMP.err && \
tophat --output-dir $TOPHAT_INDEX_DIR/tophat_out --GTF $TOPHAT_INDEX_DIR/$GTF --transcriptome-index=$TOPHAT_INDEX_PREFIX $BOWTIE2_INDEX_PREFIX > \$LOG 2> \$ERR && \
chmod -R ug+rwX,o+rX \$TOPHAT_INDEX_DIR \$LOG \$ERR"
      cmd_or_job BOWTIE2_TOPHAT_CMD 4
    else
      echo
      echo "Genome Bowtie 2 index and gtf TopHat index up to date... skipping"
      echo
    fi
  else
    if ! is_up2date $BOWTIE2_INDEX_PREFIX.[1-4].bt2 $BOWTIE2_INDEX_PREFIX.rev.[12].bt2
    then
      echo
      echo "Creating genome Bowtie 2 index..."
      echo
      BOWTIE2_CMD="\
mkdir -p $BOWTIE2_INDEX_DIR && \
ln -s -f -t $BOWTIE2_INDEX_DIR ../$GENOME_FASTA && \
module load $module_bowtie2 && \
LOG=$LOG_DIR/bowtie2_$TIMESTAMP.log && \
ERR=$LOG_DIR/bowtie2_$TIMESTAMP.err && \
bowtie2-build $BOWTIE2_INDEX_DIR/$GENOME_FASTA $BOWTIE2_INDEX_PREFIX > \$LOG 2> \$ERR && \
chmod -R ug+rwX,o+rX $BOWTIE2_INDEX_DIR \$LOG \$ERR"
      cmd_or_job BOWTIE2_CMD 4
    else
      echo
      echo "Genome Bowtie 2 index up to date... skipping"
      echo
    fi
  fi
}

create_star_index() {
  if is_genome_big
  then
    runThreadN=20
  else
    runThreadN=1
  fi
  # Since version 102 or so, Ensembl gtf files have an additional column inserted at column 3, this makes STAR genome generation crash and needs to be removed...
  if ! is_up2date $ANNOTATIONS_DIR/$GTF
  then
    cut -f -2,4- $ANNOTATIONS_DIR/$GTF > $ANNOTATIONS_DIR/$GTF
  fi
  for sjdbOverhang in 49 74 99 124 149
  do
    STAR_VERSION=$(echo $module_star | cut -f3 -d/)
    INDEX_DIR=$INSTALL_DIR/genome/star${STAR_VERSION}_index/$SOURCE$VERSION.sjdbOverhang$sjdbOverhang
    if ! is_up2date $INDEX_DIR/SAindex
    then
      echo
      echo "Creating STAR index with sjdbOverhang $sjdbOverhang..."
      echo
      STAR_CMD="\
mkdir -p $INDEX_DIR && \
module load $module_star && \
LOG=$LOG_DIR/star_${sjdbOverhang}_$TIMESTAMP.log && \
ERR=$LOG_DIR/star_${sjdbOverhang}_$TIMESTAMP.err && \
STAR --runMode genomeGenerate --genomeDir $INDEX_DIR --genomeFastaFiles $GENOME_DIR/$GENOME_FASTA --runThreadN $runThreadN --sjdbOverhang $sjdbOverhang --genomeSAindexNbases 4 --limitGenomeGenerateRAM 92798303616 --sjdbGTFfile $ANNOTATIONS_DIR/$GTF --outFileNamePrefix $INDEX_DIR/ > \$LOG 2> \$ERR && \
chmod -R ug+rwX,o+rX $INDEX_DIR \$LOG \$ERR"
      cmd_or_job STAR_CMD $runThreadN STAR_${SPECIES}_${ASSEMBLY}_${sjdbOverhang}_CMD 360G 
    else
      echo
      echo "STAR index with sjdbOverhang $sjdbOverhang up to date... skipping"
      echo
    fi
  done
}

create_ncrna_bwa_index() {
  INDEX_DIR=$ANNOTATIONS_DIR/ncrna_bwa_index
  if ! is_up2date $INDEX_DIR/$NCRNA.sa
  then
    echo
    echo "Creating ncRNA BWA index..."
    echo
    mkdir -p $INDEX_DIR
    ln -s -f -t $INDEX_DIR ../$NCRNA
    module load $module_bwa
    bwa index $INDEX_DIR/$NCRNA > $LOG_DIR/ncrna_bwa_$TIMESTAMP.log 2>&1
  else
    echo
    echo "ncRNA BWA index up to date... skipping"
    echo
  fi
}

create_rrna_bwa_index() {
  if is_up2date $ANNOTATIONS_DIR/$RRNA
  then
    INDEX_DIR=$ANNOTATIONS_DIR/rrna_bwa_index
    if ! is_up2date $INDEX_DIR/$RRNA.sa
    then
      echo
      echo "Creating rRNA BWA index..."
      echo
      mkdir -p $INDEX_DIR
      ln -s -f -t $INDEX_DIR ../$RRNA
      module load $module_bwa
      bwa index $INDEX_DIR/$RRNA > $LOG_DIR/rrna_bwa_$TIMESTAMP.log 2>&1
    else
      echo
      echo "rRNA BWA index up to date... skipping"
      echo
    fi
  fi
}

create_kallisto_index() {
  if is_up2date $ANNOTATIONS_DIR/$CDNA
  then
    INDEX_DIR=$ANNOTATIONS_DIR/cdna_kallisto_index
    if ! is_up2date $INDEX_DIR/$CDNA.idx
    then
      echo
      echo "Creating cDNA Kallisto index..."
      echo
      mkdir -p $INDEX_DIR
      ln -s -f -t $INDEX_DIR ../$CDNA
      module load $module_kallisto
      kallisto index -i $INDEX_DIR/$CDNA.idx $INDEX_DIR/$CDNA > $LOG_DIR/cdna_kallisto_$TIMESTAMP.err 2> $LOG_DIR/cdna_kallisto_$TIMESTAMP.log

    else
      echo
      echo "cDNA Kallisto index up to date... skipping"
      echo
    fi
  fi
}

create_transcripts2genes_file() {
  ANNOTATION_GTF=$ANNOTATIONS_DIR/$GTF
  if is_up2date $ANNOTATION_GTF
  then
    ANNOTATION_TX2GENES=$ANNOTATIONS_DIR/cdna_kallisto_index/${GTF/.gtf/.tx2gene}
    if ! is_up2date $ANNOTATION_TX2GENES
    then
      module load $module_R
      module load $module_mugqic_R_packages
      R --no-restore --no-save<<EOF
suppressPackageStartupMessages(library(rtracklayer))
print("Building transcripts2genes...")
gtf_file = "$ANNOTATION_GTF"

gtf = import(gtf_file, format = "gff2")
gtf=gtf[!is.na(gtf\$transcript_id),]
if ("transcript_version" %in% names(gtf@elementMetadata)){
    tx_id=paste(gtf\$transcript_id, gtf\$transcript_version, sep=".")
} else {
    tx_id=gtf\$transcript_id
}
tx2gene = cbind(tx_id=tx_id, gene_id=gtf\$gene_id)
tx2gene = unique(tx2gene)
tx2gene = as.data.frame(tx2gene)

write.table(x=tx2gene, file="$ANNOTATION_TX2GENES", sep="\t", col.names=T, row.names=F, quote=F)
EOF
    else
      echo
      echo "transcripts2genes file up to date... skipping"
      echo
    fi
  fi
}

create_gene_annotations() {
  ANNOTATION_PREFIX=$ANNOTATIONS_DIR/${GTF/.gtf}

  if ! is_up2date $ANNOTATION_PREFIX.genes.length.tsv $ANNOTATION_PREFIX.geneid2Symbol.tsv $ANNOTATION_PREFIX.genes.tsv
  then
    echo
    echo "Creating gene ID, symbol, length annotations from GTF..."
    echo
    cd $ANNOTATIONS_DIR
    module load $module_R
    module load $module_mugqic_R_packages
    module load $module_perl
    R --no-restore --no-save<<EOF
suppressPackageStartupMessages(library(gqSeqUtils))
gtf.fn     = "$GTF"
annotation = "$ANNOTATION_PREFIX"

genes = extractGeneAnnotationFromGtf(gtf.fn, feature.types = c("exon", "CDS"))
# NOTE: feature.types passed down to gqSeqUtils::calculateGeneLengthsFromGtf

write.table(genes[,c("featureID","gene.length"),],
  file = paste0(annotation, ".genes.length.tsv"), sep='\t', quote=F, row.names=F, col.names=F) # NOTE: this should not be necessary, all is in genes anyway

write.table(data.frame("ensembl"=genes[["featureID"]], "geneSymbol"=genes[["SYMBOL"]]),
  file = paste0(annotation, ".geneid2Symbol.tsv"), sep='\t', quote=F, row.names=F, col.names=F) # NOTE: this should not be necessary, all is in genes anyway

write.table(genes, file = paste0(annotation, ".genes.tsv"), sep='\t', quote=F, row.names=F, col.names=T)

EOF
  else
    echo
    echo "Gene ID, symbol, length annotations from GTF up to date... skipping"
    echo
  fi
}

create_gene_annotations_flat() {
  ANNOTATION_PREFIX=$ANNOTATIONS_DIR/${GTF/.gtf}

  if ! is_up2date $ANNOTATION_PREFIX.ref_flat.tsv
  then
    echo
    echo "Creating gene refFlat file from GTF..."
    echo
    cd $ANNOTATIONS_DIR
    module load $module_ucsc
    gtfToGenePred -genePredExt -geneNameAsName2 ${ANNOTATION_PREFIX}.gtf ${ANNOTATION_PREFIX}.refFlat.tmp.txt
    cut -f 12 ${ANNOTATION_PREFIX}.refFlat.tmp.txt > ${ANNOTATION_PREFIX}.refFlat.tmp.2.txt
    cut -f 1-10 ${ANNOTATION_PREFIX}.refFlat.tmp.txt > ${ANNOTATION_PREFIX}.refFlat.tmp.3.txt
    paste ${ANNOTATION_PREFIX}.refFlat.tmp.2.txt ${ANNOTATION_PREFIX}.refFlat.tmp.3.txt > ${ANNOTATION_PREFIX}.ref_flat.tsv
    rm ${ANNOTATION_PREFIX}.refFlat.tmp.txt ${ANNOTATION_PREFIX}.refFlat.tmp.2.txt ${ANNOTATION_PREFIX}.refFlat.tmp.3.txt
  else
    echo
    echo "Creating gene refFlat file from GTF up to date... skipping"
    echo
  fi
}


create_go_annotations() {

  if [ ! -z "${BIOMART_HOST:-}" ]
  then
    if ! is_up2date $ANNOTATIONS_DIR/$GO
    then
      echo
      echo "Retrieving GO annotations using BioMart..."
      echo
      module load $module_R
      module load $module_mugqic_R_packages
      module load $module_perl
      R --no-restore --no-save<<EOF
suppressPackageStartupMessages(library(gqSeqUtils))
httr::set_config(httr::config(cainfo = '/etc/ssl/certs/ca-bundle.crt'))
res = getBMsimple("$BIOMART_HOST","$BIOMART_MART","$BIOMART_DATASET",c("$BIOMART_GENE_ID", "$BIOMART_GO_ID", "$BIOMART_GO_NAME", "$BIOMART_GO_DEFINITION"))
res = res[ res[["$BIOMART_GENE_ID"]] != '' & res[["$BIOMART_GO_ID"]] != '', ]
write.table(res, file="$ANNOTATIONS_DIR/$GO", quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
EOF
    else
      echo
      echo "GO annotations up to date... skipping"
      echo
    fi
  fi
}

get_vcf_dbsnp() {
  if [[ ! -z "${VCF_URL:-}" && ! -z "${VCF_TBI_URL:-}" ]]
  then
    if ! is_up2date $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP[0-9]*.vcf.gz $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP[0-9]*.vcf.gz.tbi
    then
      echo
      echo "Retrieving dbSNP latest version from VCF..."
      echo
      # Set +e since zgrep exit code != 0 if not found
      set +e
      DBSNP_VERSION=`zcat $ANNOTATIONS_DIR/$VCF | grep -v "^#" | cut -f8 | grep -Po "dbSNP_\d+" | sed s/dbSNP_// | sort -ug | tail -1`
      set -e
      echo "Found VCF dbSNP version: $DBSNP_VERSION"
      if [ $DBSNP_VERSION ]
      then
        DBSNP_VCF=$ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP$DBSNP_VERSION.vcf.gz
        # Add dbSNP version to VCF filename
        # Remove variation sequences containing "." (columns 4 and 5 from uncommented lines) which make GATK crash
        module load $module_tabix
        zgrep -Pv "^[^#]\S*\t\S+\t\S+\t(\S*\.\S*|\S+\t\S*\.\S*)" $ANNOTATIONS_DIR/$VCF | bgzip > $DBSNP_VCF
        # Create tabix index for new VCF without "." variation sequences
        tabix -p vcf $DBSNP_VCF
      fi
    else
      echo
      echo "dbSNP VCF exists... assuming up to date and skipping"
      echo
    fi
  fi
}

copy_files() {
  # Uncompress files
  if ! is_up2date $GENOME_DIR/$GENOME_FASTA
  then
    if [[ $SOURCE == "UCSC" ]]
    then
      # Use UCSC utility program twoBitToFa to uncompress genome
      module load mugqic/ucsc
      twoBitToFa `download_path $GENOME_URL` $GENOME_DIR/$GENOME_FASTA
    else
      gunzip -c `download_path $GENOME_URL` > $GENOME_DIR/$GENOME_FASTA
    fi
  fi

  # Annotations are not installed for UCSC genomes
  if [[ $SOURCE != "UCSC" ]]
  then
    if ! is_up2date $ANNOTATIONS_DIR/$GTF ; then gunzip -c `download_path $GTF_URL` > $ANNOTATIONS_DIR/$GTF ; fi
    TRANSCRIPT_ID_GTF=$ANNOTATIONS_DIR/${GTF/.gtf/.transcript_id.gtf}
    if ! is_up2date $TRANSCRIPT_ID_GTF ; then grep -P "(^#|transcript_id)" $ANNOTATIONS_DIR/$GTF > $TRANSCRIPT_ID_GTF ; fi
    if ! is_up2date $ANNOTATIONS_DIR/$NCRNA ; then gunzip -c `download_path $NCRNA_URL` > $ANNOTATIONS_DIR/$NCRNA ; fi
    if ! is_up2date $ANNOTATIONS_DIR/$CDNA ; then gunzip -c `download_path $CDNA_URL` > $ANNOTATIONS_DIR/$CDNA ; fi

    # Create rRNA FASTA as subset of ncRNA FASTA keeping only sequences with "rRNA" (ignore case) in their descriptions
    if [ $(grep -q -i "rRNA" $ANNOTATIONS_DIR/$NCRNA)$? == 0 ]
    then
      if ! is_up2date $ANNOTATIONS_DIR/$RRNA
      then
        $grep_cvmfs -Pzoi "^>.*rRNA[^>]*" $ANNOTATIONS_DIR/$NCRNA | grep -v "^$" > $ANNOTATIONS_DIR/$RRNA
      fi
    fi

    if [ ! -z "${VCF_URL:-}" ]
    then
      if ! is_up2date $ANNOTATIONS_DIR/$VCF
      then
        cp `download_path $VCF_URL` $ANNOTATIONS_DIR/$VCF
      fi
    fi

    if [ ! -z "${VCF_TBI_URL:-}" ]
    then
      if ! is_up2date $ANNOTATIONS_DIR/$VCF.tbi
      then
        cp `download_path $VCF_TBI_URL` $ANNOTATIONS_DIR/$VCF.tbi
      fi
    fi

    get_vcf_dbsnp
  else
    if ! is_up2date $ANNOTATIONS_DIR/$GTF
    then
      echo "Could not find $ANNOTATIONS_DIR/$GTF...\nyou might consider to manually download a gtf file from UCSC table browser (http://genome.ucsc.edu/cgi-bin/hgTables)"
    else
      TRANSCRIPT_ID_GTF=$ANNOTATIONS_DIR/${GTF/.gtf/.transcript_id.gtf}
      if ! is_up2date $TRANSCRIPT_ID_GTF ; then grep -P "(^#|transcript_id)" $ANNOTATIONS_DIR/$GTF > $TRANSCRIPT_ID_GTF ; fi
    fi

    if ! is_up2date $ANNOTATIONS_DIR/$NCRNA
    then
      echo "Could not find $ANNOTATIONS_DIR/$NCRNA...\nyou might consider to use the ncrna.fa file from Ensembl... "
    else
      if [ $(grep -q -i "rRNA" $ANNOTATIONS_DIR/$NCRNA)$? == 0 ]
      then
        if ! is_up2date $ANNOTATIONS_DIR/$RRNA
        then
          $grep_cvmfs -Pzoi "^>.*rRNA[^>]*" $ANNOTATIONS_DIR/$NCRNA | grep -v "^$" > $ANNOTATIONS_DIR/$RRNA
        fi
      fi
    fi

    if ! is_up2date $ANNOTATIONS_DIR/$CDNA
    then
      echo "Could not find $ANNOTATIONS_DIR/$CDNA...\nyou might consider to use the cdna.fa file from Ensembl... "
    fi

  fi
}

create_genome_digest() {

  GENOME_DIGEST=$GENOME_DIR/genome_digest/

  declare -A enzymes=( ["DpnII"]="^GATC" ["MboI"]="^GATC" ["HindIII"]="A^AGCTT" ["NcoI"]="C^CATGG" ["Arima"]="^GATC,DpnII:G^ANTC")

  for enzyme in "${!enzymes[@]}"; do
    #echo "$enzyme - ${enzymes[$enzyme]}";
    ## hicup only accepts alphanumeric and underscores
    GENOME_DIGEST_FILE=HiCUP_Digest_${SPECIES}_${ASSEMBLY}_${enzyme}.txt

    if ! is_up2date $GENOME_DIGEST/$GENOME_DIGEST_FILE
    then
      echo
      echo "Creating ${enzyme} genome digest..."
      echo
      Digest_CMD="\
mkdir -p $GENOME_DIGEST && \
cd $GENOME_DIGEST  && \
ln -s -f -t $GENOME_DIGEST ../$GENOME_FASTA && \
module load $module_hicup $module_perl && \
LOG=$LOG_DIR/${enzyme}_digest_$TIMESTAMP.log && \
ERR=$LOG_DIR/${enzyme}_digest_$TIMESTAMP.err && \
ASSEMBLY_FOR_CMD=$(echo $ASSEMBLY | sed -e 's/[-.]/_/g')
hicup_digester --genome \$ASSEMBLY_FOR_CMD --re1 ${enzymes[$enzyme]},${enzyme} $GENOME_DIGEST/$GENOME_FASTA > \$LOG 2> \$ERR && \
mv Digest_\$ASSEMBLY_FOR_CMD_*${enzyme}_None_*.txt $GENOME_DIGEST_FILE" # && \
#chmod -R ug+rwX,o+rX $GENOME_DIGEST \$LOG \$ERR"
      cmd_or_job Digest_CMD 2
    else
      echo
      echo "${enzyme} genome digest is up to date... skipping"
      echo
    fi
  done
}

build_files() {
  # Create indexes
  create_picard_index
  create_samtools_index
  create_bwa_index
  create_genome_digest
  create_bismark_genome_reference
  create_bowtie_index
  create_bowtie2_tophat_index

  if is_up2date $ANNOTATIONS_DIR/$NCRNA
  then
    create_ncrna_bwa_index
    create_rrna_bwa_index
  else
    echo "Could not find $ANNOTATIONS_DIR/$NCRNA..."
    echo "No ncRNA bwa index will be created... this step is skipped..."
    echo "No rRNA bwa index will be created... this step is skipped..."
    echo "You might consider to use the ncrna.fa file from Ensembl, if available... "
  fi

  if is_up2date $ANNOTATIONS_DIR/$CDNA
  then
    create_kallisto_index
    create_transcripts2genes_file
  else
    echo "Could not find $ANNOTATIONS_DIR/$CDNA..."
    echo "No cDNA kallisto index will be created... this step is skipped..."
    echo "No tx2gene file will be created...this step is skipped..."
    echo "You might consider to use the cdna.fa file from Ensembl, if available... "
  fi

  if is_up2date $ANNOTATIONS_DIR/$GTF
  then
    create_star_index
    create_gene_annotations
    create_gene_annotations_flat
  else
    echo "Could not find $ANNOTATIONS_DIR/$GTF..."
    echo "No tophat index will be created... this step is skipped..."
    echo "No gene annotations will be performed... this step is skipped..."
    echo "You might consider to manually download a gtf file from UCSC table browser (http://genome.ucsc.edu/cgi-bin/hgTables)"
  fi

  # Annotations are not installed for UCSC genomes
  if [[ $SOURCE != "UCSC" ]]
  then
    create_go_annotations
  fi
}

create_genome_ini_file() {
  INI=$INSTALL_DIR/$SPECIES.$ASSEMBLY.${SOURCE}${VERSION}.ini
  echo "\
[DEFAULT]
scientific_name=$SPECIES
common_name=$COMMON_NAME
assembly=$ASSEMBLY
assembly_synonyms=$ASSEMBLY_SYNONYMS
source=$SOURCE
version=$VERSION
module_star=$module_star" > $INI

  if [ ! -z "${DBSNP_VERSION:-}" ]
  then
    mv $INI ${INI/ini/dbSNP${DBSNP_VERSION}.ini}
    INI=${INI/ini/dbSNP${DBSNP_VERSION}.ini}
    echo "\
dbsnp_version=$DBSNP_VERSION" >> $INI
  fi
}

install_genome() {

  SPECIES="$1"
  COMMON_NAME="$2"
  ASSEMBLY="$3"
  ASSEMBLY_SYNONYMS="$4"
  SOURCE="$5"
  VERSION="$6"

  init_install
  set_urls
  download_urls
  copy_files
  build_files
  create_genome_ini_file

  # Add permissions
  chmod -R ug+rwX,o+rX $INSTALL_DIR
}

