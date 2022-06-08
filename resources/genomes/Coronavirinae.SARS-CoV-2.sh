#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Coronavirinae
COMMON_NAME="Coronavirus"
ASSEMBLY=SARS-CoV2
ASSEMBLY_SYNONYMS=MN908947.3
SOURCE=GenBank
VERSION=2020-02-17
BIOMART_HOST=

module_snpeff=mugqic/snpEff/4.3
module_tabix=mugqic/tabix/0.2.6
module_java=mugqic/java/openjdk-jdk1.8.0_72

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

# Overwrite install_genome since GenBank/NCBI genome is used instead of Ensembl
install_genome() {

  SPECIES="$1"
  COMMON_NAME="$2"
  ASSEMBLY="$3"
  ASSEMBLY_SYNONYMS="$4"
  SOURCE="$5"
  VERSION="$6"

  init_install
  set_urls

  # NCBI no alternate genome is used since Ensembl version contains all alternate haplotypes, does not contain EBV, and chromosome entries are not sorted (which causes problem for GATK)
  GENOME_URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz

  download_urls
  # set +e since gunzip human_g1k_v37.fasta.gz exit code != 0 ("gzip: human_g1k_v37.fasta.gz: decompression OK, trailing garbage ignored")
  set +e
  copy_files
  build_files
  create_genome_ini_file

  # Add permissions
  chmod -R ug+rwX,o+rX $INSTALL_DIR
}

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
