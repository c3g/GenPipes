#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Macaca_mulatta
COMMON_NAME="Macaque"
ASSEMBLY=MMUL_1
ASSEMBLY_SYNONYMS=rheMac3
SOURCE=Ensembl
VERSION=84
BIOMART_HOST=mar2016.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly

