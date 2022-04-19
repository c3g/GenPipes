#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Pocillopora_damicornis
COMMON_NAME="Lace Coral"
ASSEMBLY=GCF_003704095.1
ASSEMBLY_SYNONYMS=
SOURCE=UCSC
VERSION=2019-08-07

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

set_urls() {
    GENOME_URL=https://hgdownload.soe.ucsc.edu/hubs/GCF/003/704/095/$ASSEMBLY/$ASSEMBLY.2bit
}

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
