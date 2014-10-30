#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Homo_sapiens
COMMON_NAME="Human"
ASSEMBLY=GRCh37
ASSEMBLY_SYNONYMS=hg19
SOURCE=Ensembl
VERSION=75

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

# Download dbSNP directly from NCBI since it is more up to date
get_vcf_dbsnp() {
  DBSNP_VERSION=142
  DBSNP_URL=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b${DBSNP_VERSION}_GRCh37p13/VCF/All.vcf.gz
  download_url $DBSNP_URL
  download_url $DBSNP_URL.tbi
  cp `download_path $DBSNP_URL` $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP$DBSNP_VERSION.vcf.gz
  cp `download_path $DBSNP_URL.tbi` $ANNOTATIONS_DIR/$SPECIES.$ASSEMBLY.dbSNP$DBSNP_VERSION.vcf.gz.tbi
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

  # 1000Genomes genome is used since Ensembl version chromosome entries are not sorted (which causes problem for GATK)
  GENOME_URL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

  download_urls
  # set +e since gunzip human_g1k_v37.fasta.gz exit code != 0 ("gzip: human_g1k_v37.fasta.gz: decompression OK, trailing garbage ignored")
  set +e
  copy_files
  set -e
  build_files
  create_genome_ini_file

  # Add permissions
  chmod -R ug+rwX,o+rX $INSTALL_DIR
}

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
