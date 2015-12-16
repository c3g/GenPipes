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
source $GENOME_INSTALL_SCRIPT_DIR/../install_genome.sh

install_cosmic() {
  # download COSMIC manually using your username pass
  # and put them in DOWNLOAD_DIR
  # cd $DOWNLOAD_DIR
  # lftp -e "get CosmicNonCodingVariants.vcf.gz -o CosmicNonCodingVariants.72.vcf.gz" -u "<USER>,<PASS>" sftp://sftp-cancer.sanger.ac.uk/files/grch37/cosmic/v72/VCF/
  # lftp -e "get CosmicCodingMuts.vcf.gz -o CosmicCodingMuts.72.vcf.gz" -u "<USER>,<PASS>" sftp://sftp-cancer.sanger.ac.uk/files/grch37/cosmic/v72/VCF/

  module load mugqic/htslib/1.2.1
  cat <(zgrep -v "^#" $DOWNLOAD_DIR/CosmicCodingMuts.72.vcf.gz) <(zgrep -v "^#" $DOWNLOAD_DIR/CosmicNonCodingVariants.72.vcf.gz) | sort -S2G -k1V,1V -k2n,2n | cat <(zgrep "^#" $DOWNLOAD_DIR/CosmicCodingMuts.72.vcf.gz) - | bgzip -c > $ANNOTATIONS_DIR/cosmic.72.vcf.gz
  tabix -p vcf $ANNOTATIONS_DIR/cosmic.72.vcf.gz
}

install_genome() {

  SPECIES="$1"
  COMMON_NAME="$2"
  ASSEMBLY="$3"
  ASSEMBLY_SYNONYMS="$4"
  SOURCE="$5"
  VERSION="$6"

  INSTALL_HOME=$MUGQIC_INSTALL_HOME_PRIVATE
  init_install

  install_cosmic

  # Add permissions
  chmod -R ug+rwX,o+rX $INSTALL_DIR
}

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
