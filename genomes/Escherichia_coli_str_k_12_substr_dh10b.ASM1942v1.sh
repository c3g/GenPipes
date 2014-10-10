#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Escherichia_coli_str_k_12_substr_dh10b
COMMON_NAME="Escherichia coli strain K-12 substrain DH10B"
ASSEMBLY=ASM1942v1
ASSEMBLY_SYNONYMS=
SOURCE=EnsemblGenomes
VERSION=23

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
