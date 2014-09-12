#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SPECIES=Homo_sapiens
COMMON_NAME="Human"
ASSEMBLY=GRCh38
ASSEMBLY_SYNONYMS=hg38
SOURCE=Ensembl
RELEASE=76

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $GENOME_INSTALL_SCRIPT_DIR/install_genome.sh

# Overwrite install_genome since NCBI genome is used instead of Ensembl
install_genome() {

  SPECIES="$1"
  COMMON_NAME="$2"
  ASSEMBLY="$3"
  ASSEMBLY_SYNONYMS="$4"
  SOURCE="$5"
  RELEASE="$6"

  init_install
  set_urls

  # NCBI no alternate genome is used since Ensembl version contains all alternate haplotypes, does not contain EBV, and chromosome entries are not sorted (which causes problem for GATK)
  GENOME_URL=ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

  download_urls
  copy_files

  echo Update Ensembl GTF to match NCBI genome...
  # Remove Ensembl GTF haplotype annotations, adjust mitochondria and "GK" annotation version names
  grep -v "^CHR_H" $ANNOTATIONS_DIR/$GTF | sed 's/^MT\t/M\t/' | perl -pe "s/^([GK]\S+)\.(\d+)\t/\1v\2\t/" > $ANNOTATIONS_DIR/$GTF.tmp
  # Update Ensembl GTF annotation IDs to match NCBI genome chromosome IDs
  grep "^>" $GENOME_DIR/$GENOME_FASTA | cut -f1 -d\  | cut -c 2- | perl -pe 's/^(chr([^_\n]*))$/\1\t\2/' | perl -pe 's/^(chr[^_]*_([^_\n]*)(_\S+)?)$/\1\t\2/' | awk -F"\t" 'FNR==NR{id[$2]=$1; next}{OFS="\t"; if (id[$1]) {print id[$1],$0} else {print $0}}' - $ANNOTATIONS_DIR/$GTF.tmp | cut -f1,3- > $ANNOTATIONS_DIR/$GTF
  rm $ANNOTATIONS_DIR/$GTF.tmp

  build_files
  create_genome_ini_file

  # Add permissions
  chmod ug+rwX,o+rX $INSTALL_DIR
}

install_genome "$SPECIES" "$COMMON_NAME" "$ASSEMBLY" "$ASSEMBLY_SYNONYMS" "$SOURCE" "$RELEASE"

################################################################################
# Write below all commands to install additional data files specific to this genome assembly
