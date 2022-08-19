#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Anas_platyrhynchos_platyrhynchos
COMMON_NAME="Duck"
ASSEMBLY=CAU_duck1.0
ASSEMBLY_SYNONYMS=
SOURCE=Ensembl
VERSION=98
BIOMART_HOST=sep2019.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
