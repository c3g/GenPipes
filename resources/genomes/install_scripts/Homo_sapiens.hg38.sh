#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Homo_sapiens
COMMON_NAME="Human"
ASSEMBLY=hg38
ASSEMBLY_SYNONYMS=GRCh38
SOURCE=UCSC
VERSION=2014-01-09

module_snpeff=mugqic/snpEff/4.2
module_tabix=mugqic/tabix/0.2.6
module_java=mugqic/java/openjdk-jdk1.8.0_72

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

# Download dbSNP directly from NCBI since it is more up to date
get_vcf_dbsnp() {
  DBSNP_VERSION=142
  DBSNP_URL=ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/archive/human_9606_b${DBSNP_VERSION}_GRCh38/VCF/All.vcf.gz
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
    DBNSFP_URL=ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.1c.zip
    DBSNSFP_VERSION=dbNSFPv3.1c
    DBSNSFP=$ANNOTATIONS_DIR/$DBSNSFP_VERSION/$DBSNSFP_VERSION.txt.gz
    if ! is_up2date $DBSNSFP.txt.gz
    then
      mkdir -p $ANNOTATIONS_DIR/$DBSNSFP_VERSION/
      if ! is_up2date `download_path $DBNSFP_URL`; then
        download_url $DBNSFP_URL
        cp dbnsfp.softgenetics.com/dbNSFPv3.1c.zip $ANNOTATIONS_DIR/$DBSNSFP_VERSION/
      fi
      unzip $ANNOTATIONS_DIR/$DBSNSFP_VERSION/$DBSNSFP_VERSION.zip -d $ANNOTATIONS_DIR/$DBSNSFP_VERSION/
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
      java -Xmx8G -jar $SNPEFF_HOME/SnpSift.jar dbnsfp -v -db $DBSNSFP.txt.gz $DBSNP > $DBSNP_ANNOTATED.tmp
      # Filter autosomal SNPs
      (grep "^#" $DBSNP_ANNOTATED.tmp > $DBSNP_ANNOTATED; cat $DBSNP_ANNOTATED.tmp | grep -v "^#" | awk 'BEGIN{FS="\t";OFS="\t";}{ if ( $1 ~ /^[0-9]/) {print $0;}}' >> $DBSNP_ANNOTATED )
      for POP_FREQ in 1000Gp1_EUR_AF 1000Gp1_AFR_AF 1000Gp1_ASN_AF;
      do
        cat $DBSNP_ANNOTATED | sed -e 's/dbNSFP_'$POP_FREQ'/AF/g' > $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP${DBSNP_VERSION}_${POP_FREQ}.vcf
        module load $module_tabix
        bgzip $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP${DBSNP_VERSION}_${POP_FREQ}.vcf
        tabix -p vcf $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP${DBSNP_VERSION}_${POP_FREQ}.vcf.gz
      done
    fi
    # set the default allele frequency for a population (hapmap CEU)
    population_AF=1000Gp1_EUR_AF
}

# Overwrite install_genome since NCBI genome is used instead of Ensembl
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
