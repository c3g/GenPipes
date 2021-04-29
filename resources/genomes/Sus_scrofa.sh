#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Sus_scrofa
COMMON_NAME="Pig"
ASSEMBLY=Sscrofa11.1
ASSEMBLY_SYNONYMS=susScr11
SOURCE=Ensembl
VERSION=99
BIOMART_HOST=jan2020.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh $@

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
