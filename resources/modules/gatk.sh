#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=GenomeAnalysisTK

# What follows stands for version <= 3.8
#VERSION=3.8
#ARCHIVE=$SOFTWARE-$VERSION.tar.bz2
#echo "Prior to install the gatk module, you must download the archive $ARCHIVE manually, if not done already, from http://www.broadinstitute.org/gatk/download since it requires a license agreement.
#Once downloaded, copy it in \$MUGQIC_INSTALL_HOME_DEV/archive/ or \$MUGQIC_INSTALL_HOME/archive/"
#SOFTWARE_DIR=$SOFTWARE-$VERSION
#ARCHIVE_URL=https://www.broadinstitute.org/gatk/download/auth?package=GATK
#JAR=$SOFTWARE.jar

# What follows stands for version > 3.8
VERSION=4.1.0.0
ARCHIVE=gatk-${VERSION}.zip
SOFTWARE_DIR=${SOFTWARE}-${VERSION}
ARCHIVE_URL=https://github.com/broadinstitute/gatk/releases/download/${VERSION}/gatk-${VERSION}.zip
JAR=gatk-package-${VERSION}-local.jar


build() {
# What follows stands for version <= 3.8
#  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
#  tar -jxvf $INSTALL_DOWNLOAD/$ARCHIVE --directory=$INSTALL_DIR/$SOFTWARE_DIR

# What follows stands for version > 3.8
  mkdir -p $INSTALL_DIR
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE -d $INSTALL_DIR/$SOFTWARE_DIR
  mv $INSTALL_DIR/$SOFTWARE_DIR/${ARCHIVE//.zip/}/* $INSTALL_DIR/$SOFTWARE_DIR/
  rm -rf $INSTALL_DIR/$SOFTWARE_DIR/${ARCHIVE//.zip/}
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH	            \$root
setenv          GATK_JAR            \$root/$JAR
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
