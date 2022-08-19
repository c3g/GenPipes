#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Homo_sapiens
COMMON_NAME="Human"
ASSEMBLY=hg19
ASSEMBLY_SYNONYMS=GRCh37
SOURCE=UCSC
VERSION=2009-03-08
DBSNP_VERSION=150
ONEKG_PHASE="p3"

module_snpeff=mugqic/snpEff/4.2
module_tabix=mugqic/tabix/0.2.6
module_java=mugqic/java/openjdk-jdk1.8.0_72

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

# Download dbSNP directly from NCBI since it is more up to date
get_vcf_dbsnp() {
  DBSNP_VERSION=150
#  DBSNP_URL=ftp://ftp.ncbi.nih.gov/snp/organisms/archive/human_9606_b${DBSNP_VERSION}_GRCh37p13/VCF/GATK/00-All.vcf.gz    # to use with DBSNP_VERSION = 142 & 144
  DBSNP_URL=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b${DBSNP_VERSION}_GRCh37p13/VCF/GATK/00-All.vcf.gz            # to use with DBSNP_VERSION >= 146
  DBSNP=$ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP$DBSNP_VERSION.vcf.gz

  if ! is_up2date $DBSNP $DBSNP.tbi
  then
    download_url $DBSNP_URL
    download_url $DBSNP_URL.tbi
    cp `download_path $DBSNP_URL` $DBSNP
    cp `download_path $DBSNP_URL.tbi` $DBSNP.tbi
  else
    echo
    echo "dbSNP file $DBSNP up to date... skipping"
  fi
}

# Download dbNSFP and generate vcfs required to run VerifyBamId
get_dbNSFP() {
    DBSNSFP_VERSION=dbNSFPv3.5a
    DBNSFP_URL=ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/${DBSNSFP_VERSION}.zip
    DBSNSFP=$ANNOTATIONS_DIR/$DBSNSFP_VERSION/$DBSNSFP_VERSION
    if ! is_up2date $DBSNSFP.txt.gz
        then
        mkdir -p $ANNOTATIONS_DIR/$DBSNSFP_VERSION/
        if ! is_up2date `download_path $DBNSFP_URL`; then
            download_url $DBNSFP_URL
            cp dbnsfp.softgenetics.com/${DBSNSFP_VERSION}.zip $ANNOTATIONS_DIR/$DBSNSFP_VERSION/
        fi
        unzip $ANNOTATIONS_DIR/$DBSNSFP_VERSION/$DBSNSFP_VERSION.zip -d $ANNOTATIONS_DIR/$DBSNSFP_VERSION/
        rm $ANNOTATIONS_DIR/$DBSNSFP_VERSION/$DBSNSFP_VERSION.zip
        (head -n 1 $ANNOTATIONS_DIR/$DBSNSFP_VERSION/*_variant.chr1 ; cat $ANNOTATIONS_DIR/$DBSNSFP_VERSION/*_variant.chr* | grep -v "^#" ) > $DBSNSFP.txt
        module load $module_tabix
        bgzip $DBSNSFP.txt
        tabix -s 1 -b 2 -e 2 $DBSNSFP.txt.gz
        rm $ANNOTATIONS_DIR/$DBSNSFP_VERSION/*_variant.chr*
    fi
    # Extract allelic frequencies for HAPMAP human populations and annotate dbsnp VCF
    DBSNP_ANNOTATED=$ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP${DBSNP_VERSION}_annotated.vcf
    if ! is_up2date $DBSNP_ANNOTATED; then
        module load $module_snpeff $module_java
        java -Xmx8G -jar $SNPEFF_HOME/SnpSift.jar dbnsfp -v -db $DBSNSFP.txt.gz $DBSNP > $DBSNP_ANNOTATED
        for POP_FREQ in 1000G${ONEKG_PHASE}_EUR_AF 1000G${ONEKG_PHASE}_AFR_AF 1000G${ONEKG_PHASE}_ASN_AF;
        do
            cat $DBSNP_ANNOTATED | sed -e 's/dbNSFP_'$POP_FREQ'/AF/g' > $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP${DBSNP_VERSION}_${POP_FREQ}.vcf
            bgzip $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP${DBSNP_VERSION}_${POP_FREQ}.vcf
            tabix -s 1 -b 2 -e 2 $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP${DBSNP_VERSION}_${POP_FREQ}.vcf.gz
        done
    else
        echo "$DBSNP_ANNOTATED is up to date... skipping"
    fi
    # set the default allele frequency for a population (hapmap CEU)
    population_AF=1000Gp1_EUR_AF
}

# rewrite install_genome to add dbsnp/dbnsfp annotations
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
  # set +e since gunzip human_g1k_v37.fasta.gz exit code != 0 ("gzip: human_g1k_v37.fasta.gz: decompression OK, trailing garbage ignored")
  set +e
  copy_files
  get_vcf_dbsnp  
  get_dbNSFP
  set -e
  build_files
  create_genome_ini_file

  # Add permissions
  chmod -R ug+rwX,o+rX $INSTALL_DIR
}

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
