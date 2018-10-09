#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Arabidopsis_thaliana
COMMON_NAME="Arabidopsis thaliana,Mouse-ear cress,Thale cress"
ASSEMBLY=TAIR10
ASSEMBLY_SYNONYMS=
SOURCE=EnsemblGenomes
VERSION=23
BIOMART_HOST=plants.ensembl.org
BIOMART_MART=plants_mart

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
