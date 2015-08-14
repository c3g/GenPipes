#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

GENOME_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

$GENOME_INSTALL_SCRIPT_DIR/Escherichia_coli_str_k_12_substr_dh10b.ASM1942v1.sh
$GENOME_INSTALL_SCRIPT_DIR/Saccharomyces_cerevisiae.R64-1-1.sh
$GENOME_INSTALL_SCRIPT_DIR/Schizosaccharomyces_pombe.ASM294v2.sh
$GENOME_INSTALL_SCRIPT_DIR/Caenorhabditis_elegans.WBcel235.sh
$GENOME_INSTALL_SCRIPT_DIR/Arabidopsis_thaliana.TAIR10.sh
$GENOME_INSTALL_SCRIPT_DIR/Drosophila_melanogaster.BDGP5.sh
$GENOME_INSTALL_SCRIPT_DIR/Gallus_gallus.Galgal4.sh
$GENOME_INSTALL_SCRIPT_DIR/Canis_familiaris.CanFam3.1.sh
$GENOME_INSTALL_SCRIPT_DIR/Bos_taurus.UMD3.1.sh
$GENOME_INSTALL_SCRIPT_DIR/Mus_musculus.NCBIM37.sh
$GENOME_INSTALL_SCRIPT_DIR/Mus_musculus.GRCm38.sh
$GENOME_INSTALL_SCRIPT_DIR/Mus_musculus.mm9.sh
$GENOME_INSTALL_SCRIPT_DIR/Mus_musculus.mm10.sh
$GENOME_INSTALL_SCRIPT_DIR/Rattus_norvegicus.Rnor_5.0.sh
$GENOME_INSTALL_SCRIPT_DIR/Rattus_norvegicus.rn5.sh
$GENOME_INSTALL_SCRIPT_DIR/Homo_sapiens.GRCh38.sh
$GENOME_INSTALL_SCRIPT_DIR/Homo_sapiens.GRCh37.sh
$GENOME_INSTALL_SCRIPT_DIR/Homo_sapiens.hg19.sh
