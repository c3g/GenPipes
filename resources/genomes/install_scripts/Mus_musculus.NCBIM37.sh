#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Mus_musculus
COMMON_NAME="Mouse"
ASSEMBLY=NCBIM37
ASSEMBLY_SYNONYMS=mm9
SOURCE=Ensembl
VERSION=67
BIOMART_HOST=may2012.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
