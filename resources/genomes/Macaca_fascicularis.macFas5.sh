#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Macaca_fascicularis
COMMON_NAME="Macaque"
ASSEMBLY=macFas5
ASSEMBLY_SYNONYMS=Macaca_fascicularis_5.0
SOURCE=UCSC
VERSION=2013-06-27

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly

