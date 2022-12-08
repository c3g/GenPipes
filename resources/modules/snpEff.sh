#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=snpEff
VERSION=5.1
# Replace "." in official version number by "_" in archive version number
ARCHIVE=${SOFTWARE}_v${VERSION//./_}_core.zip
ARCHIVE_URL=https://snpeff.blob.core.windows.net/versions/${ARCHIVE}
SOFTWARE_DIR=${SOFTWARE}_${VERSION//./_}

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE

  # Install databases
  echo "Installing GRCh37.87 database"
#  java -jar $SOFTWARE/snpEff.jar download GRCh37.75 -verbose
  java -jar $SOFTWARE/snpEff.jar download GRCh37.87 -verbose
  echo "Installing GRCh38.99 database"
#  java -jar $SOFTWARE/snpEff.jar download GRCh38.86 -verbose
  java -jar $SOFTWARE/snpEff.jar download GRCh38.99 -verbose
  echo "Installing GRCm38.99 database"
#  java -jar $SOFTWARE/snpEff.jar download GRCm38.86 -verbose  
  java -jar $SOFTWARE/snpEff.jar download GRCm38.99 -verbose
  echo "Installing hg19 database"
  java -jar $SOFTWARE/snpEff.jar download hg19 -verbose
  echo "Installing hg38 database"
  java -jar $SOFTWARE/snpEff.jar download hg38 -verbose
  echo "Installing mm10 database"
  java -jar $SOFTWARE/snpEff.jar download mm10 -verbose

  # Install software
  mv -i $SOFTWARE $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          SNPEFF_HOME         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
