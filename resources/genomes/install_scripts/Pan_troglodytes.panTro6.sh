#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Pan_troglodytes
COMMON_NAME="Chimpanzee"
ASSEMBLY=panTro6
ASSEMBLY_SYNONYMS=Clint_PTRv2
SOURCE=UCSC
VERSION=2018-03-24

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly

