#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Homo_sapiens
COMMON_NAME="Human"
ASSEMBLY=GRCh37
ASSEMBLY_SYNONYMS=hg19
SOURCE=Ensembl
RELEASE=75

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

# Overwrite install_genome since NCBI genome is used instead of Ensembl
install_genome() {

  SPECIES="$1"
  COMMON_NAME="$2"
  ASSEMBLY="$3"
  ASSEMBLY_SYNONYMS="$4"
  SOURCE="$5"
  RELEASE="$6"

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
  chmod ug+rwX,o+rX $INSTALL_DIR
}

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$RELEASE"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
