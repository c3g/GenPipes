#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Pongo_abelii
COMMON_NAME="Orangutan"
ASSEMBLY=ponAbe3
ASSEMBLY_SYNONYMS=Susie_PABv2
SOURCE=UCSC
VERSION=2018-03-26

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly

