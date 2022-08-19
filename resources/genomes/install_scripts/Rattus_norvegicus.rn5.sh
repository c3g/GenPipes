#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Rattus_norvegicus
COMMON_NAME="Rat"
ASSEMBLY=rn5
ASSEMBLY_SYNONYMS=Rnor_5.0
SOURCE=UCSC
VERSION=2012-03-23

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
