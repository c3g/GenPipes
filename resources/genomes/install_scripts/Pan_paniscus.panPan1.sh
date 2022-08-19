#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Pan_paniscus
COMMON_NAME="Bonobo"
ASSEMBLY=panPan1
ASSEMBLY_SYNONYMS=panpan1.1
SOURCE=UCSC
VERSION=2013-03-07

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
