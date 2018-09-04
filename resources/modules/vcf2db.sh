#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=vcf2db
VERSION=master_20180427
ARCHIVE=$SOFTWARE-${VERSION}.zip
ARCHIVE_URL=https://github.com/quinlan-lab/vcf2db/archive/master.zip
SOFTWARE_DIR=$SOFTWARE-${VERSION}

build() {
  cd $INSTALL_DIR
  git clone https://github.com/quinlan-lab/vcf2db
  mv vcf2db $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  module load mugqic/anaconda/2-4.0.0
  conda install -y gcc snappy # install the C library for snappy
  conda install -c conda-forge python-snappy
  conda install -c bioconda cyvcf2 peddy
  module unload module load mugqic/anaconda/2-4.0.0
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
