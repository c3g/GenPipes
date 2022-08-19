#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Caenorhabditis_elegans
COMMON_NAME="Caenorhabditis elegans,Roundworm"
ASSEMBLY=WBcel235
ASSEMBLY_SYNONYMS=
SOURCE=Ensembl
VERSION=77
BIOMART_HOST=oct2014.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
