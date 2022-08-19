#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Gorilla_gorilla
COMMON_NAME="Gorilla"
ASSEMBLY=gorGor3
ASSEMBLY_SYNONYMS=gorGor3.1
SOURCE=UCSC
VERSION=2011-10-14

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
