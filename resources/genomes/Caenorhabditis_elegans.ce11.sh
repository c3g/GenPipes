#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Caenorhabditis_elegans
COMMON_NAME="Nematodes"
ASSEMBLY=ce11
ASSEMBLY_SYNONYMS=WS220
SOURCE=UCSC
VERSION=2015-06-10
BIOMART_HOST=metazoa.ensembl.org
BIOMART_MART=metazoa_mart

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
