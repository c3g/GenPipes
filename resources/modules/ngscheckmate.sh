#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=NGSCheckMate
VERSION=1.0.0_rjme
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/robrollback/${SOFTWARE}/archive/v${VERSION%_*}.tar.gz
#ARCHIVE_URL=https://github.com/parklab/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  git clone https://github.com/robrollback/${SOFTWARE}.git 

  cat <<EOF > $SOFTWARE/ncm.conf
SAMTOOLS=samtools
BCFTOOLS=bcftools
REF=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa
EOF

  # Install software
  mv -i $SOFTWARE $INSTALL_DIR/$SOFTWARE_DIR
  chmod a+x -R $INSTALL_DIR/$SOFTWARE_DIR/*
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/
setenv          NCM_HOME            \$root/
setenv          CHECKMATE_PATH      \$root/
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
