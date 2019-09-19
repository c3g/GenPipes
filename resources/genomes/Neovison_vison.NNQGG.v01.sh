#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Neovison_vison
COMMON_NAME="American Mink"
ASSEMBLY=NNQGG.v01
ASSEMBLY_SYNONYMS=
SOURCE=Ensembl
VERSION=96
BIOMART_HOST=apr2019.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
