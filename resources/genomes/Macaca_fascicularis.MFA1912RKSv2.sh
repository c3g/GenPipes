#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Macaca_fascicularis
COMMON_NAME="Crab-eating Macaque"
ASSEMBLY=MFA1912RKSv2
ASSEMBLY_SYNONYMS=
SOURCE=NCBI
VERSION=2021-12-09
ACCESSION_ID=GCF_012559485.2 
URL_PREFIX=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/012/559/485/${ACCESSION_ID}_${ASSEMBLY}/${ACCESSION_ID}_${ASSEMBLY}

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
