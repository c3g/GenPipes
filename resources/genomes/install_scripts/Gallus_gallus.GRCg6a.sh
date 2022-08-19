#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Gallus_gallus
COMMON_NAME="Chicken"
ASSEMBLY=GRCg6a
ASSEMBLY_SYNONYMS=galGal6
SOURCE=Ensembl
VERSION=98
BIOMART_HOST=sep2019.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
