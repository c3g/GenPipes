#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Pan_troglodytes
COMMON_NAME="Chimpanzee"
ASSEMBLY=panTro4
ASSEMBLY_SYNONYMS=CHIMP2.1.4
SOURCE=UCSC
VERSION=2012-01-09

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly

