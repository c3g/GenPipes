#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a genome install script template which should be copied and used for
# consistency between genome paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


SPECIES=Apis_mellifera  ## TO BE MODIFIED WITH e.g. Homo_sapiens, Danio_rerio, Drosophila_melanogaster, etc.
COMMON_NAME="European honey bee"  ## TO BE MODIFIED WITH COMMA-SEPARATED LIST OF COMMON NAMES e.g. "Human", "Zebrafish", "Fruit fly", etc. (WITH "" IF SPACES)
ASSEMBLY=Amel_4.5  ## TO BE MODIFIED WITH e.g. GRCh37, Zv9, BDGP5, etc.
ASSEMBLY_SYNONYMS=  ## TO BE MODIFIED WITH COMMA-SEPARATED LIST OF UCSC/Ensembl OR OTHER ASSEMBLY SYNONYMS e.g. hg19, danRer7, dm3, etc. OR LEFT BLANK (DO NOT DELETE THIS LINE THOUGH)
SOURCE=EnsemblGenomes  ## TO BE MODIFIED WITH SPECIFIC SOURCE
VERSION=34  ## TO BE MODIFIED WITH SOURCE VERSION
BIOMART_HOST=metazoa.ensembl.org
BIOMART_MART=metazoa_mart_34

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" ; echo -e $GENOME_INSTALL_SCRIPT_DIR
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly