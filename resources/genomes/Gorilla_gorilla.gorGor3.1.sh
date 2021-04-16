#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Gorilla_gorilla
COMMON_NAME="Gorilla"
ASSEMBLY=gorGor3.1
ASSEMBLY_SYNONYMS=Gorgor3
SOURCE=Ensembl
VERSION=90

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
