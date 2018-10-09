#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Danio_rerio
COMMON_NAME="Zebrafish"
ASSEMBLY=GRCz11
ASSEMBLY_SYNONYMS=danRer11
SOURCE=Ensembl
VERSION=92
BIOMART_HOST=apr2018.archive.ensembl.org

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
