#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Macaca_fascicularis
COMMON_NAME="Crab-eating Macaque"
ASSEMBLY=Macaca_fascicularis_6.0
ASSEMBLY_SYNONYMS=
SOURCE=Ensembl
VERSION=106
BIOMART_HOST=apr2022.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly

