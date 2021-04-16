#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=hicpro
VERSION=2.11.1
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/nservant/HiC-Pro/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=HiC-Pro_$VERSION

build() {

  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd HiC-Pro-$VERSION

  module load mugqic/bowtie2/2.3.5 mugqic/samtools/1.9 mugqic/R_Bioconductor/3.5.3_3.8 mugqic/python/2.7.16
  sed -i "s|PREFIX = |PREFIX = ${INSTALL_DIR}|" config-install.txt
  sed -i "s|BOWTIE2_PATH = |BOWTIE2_PATH = $(dirname $(which bowtie2))|" config-install.txt
  sed -i "s|SAMTOOLS_PATH = |SAMTOOLS_PATH = $(dirname $(which samtools))|" config-install.txt
  sed -i "s|R_PATH = |R_PATH = $(dirname $(which R))|" config-install.txt
  sed -i "s|PYTHON_PATH = |PYTHON_PATH = $(dirname $(which python))|" config-install.txt
  sed -i "s|CLUSTER_SYS = |CLUSTER_SYS = SLURM|" config-install.txt

  make configure
  make install
}


#Module definition to use
module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"HiCPlotter for interaction matrix visualization for Hi-C\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
"
}


# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
