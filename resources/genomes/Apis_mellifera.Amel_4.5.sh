#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Apis_mellifera
COMMON_NAME="Honey bee"
ASSEMBLY=Amel_4.5
ASSEMBLY_SYNONYMS=
SOURCE=EnsemblGenomes
VERSION=34
BIOMART_HOST=metazoa.ensembl.org
BIOMART_MART=metazoa_mart

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" ; echo -e $GENOME_INSTALL_SCRIPT_DIR
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

set_urls() {
    RELEASE_URL=ftp://ftp.ensemblgenomes.org/pub/release-$VERSION

    # Retrieve Ensembl Genomes species information
    SPECIES_URL=$RELEASE_URL/species.txt
    download_url $SPECIES_URL
    SPECIES_LINE="`awk -F"\t" -v species=${SPECIES,,} '$2 == species' $(download_path $SPECIES_URL)`"

    # Retrieve species division (Bacteria|Fungi|Metazoa|Plants|Protists)
    DIVISION=`echo "$SPECIES_LINE" | cut -f3 | sed "s/^Ensembl//"`

    EG_SPECIES=${SPECIES,,}
    EG_BASENAME=$SPECIES.`echo "$SPECIES_LINE" | cut -f6`

    URL_PREFIX=$RELEASE_URL/${DIVISION,}
    GENOME_URL=$URL_PREFIX/fasta/$EG_SPECIES/dna/$EG_BASENAME.dna.toplevel.fa.gz
    NCRNA_URL=$URL_PREFIX/fasta/$EG_SPECIES/ncrna/$EG_BASENAME.ncrna.fa.gz
    CDNA_URL=$URL_PREFIX/fasta/$EG_SPECIES/cdna/$EG_BASENAME.cdna.all.fa.gz
    GTF_URL=$URL_PREFIX/gtf/$EG_SPECIES/$EG_BASENAME.$VERSION.gtf.gz
    if [[ `echo "$SPECIES_LINE" | cut -f8` == "Y"  ]]
    then
      VCF_URL=$URL_PREFIX/vcf/${SPECIES,,}/${SPECIES,,}.vcf.gz
    else
      echo "VCF not available for $SPECIES"
    fi
    # Retrieve species short name e.g. "athaliana" for "Arabidopsis_thaliana"
    SPECIES_SHORT_NAME=`echo ${SPECIES:0:1}${SPECIES#*_} | tr [:upper:] [:lower:]`
    BIOMART_DATASET=${SPECIES_SHORT_NAME}_eg_gene
    BIOMART_GENE_ID=ensembl_gene_id
    BIOMART_GO_ID=go_id
    BIOMART_GO_NAME=name_1006
    BIOMART_GO_DEFINITION=definition_1006
}


install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$VERSION"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
