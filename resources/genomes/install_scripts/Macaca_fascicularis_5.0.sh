#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Macaca_fascicularis
COMMON_NAME="Macaque"
ASSEMBLY=Macaca_fascicularis_5.0
ASSEMBLY_SYNONYMS=macFas5
SOURCE=Ensembl
VERSION=102
BIOMART_HOST=nov2020.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly

