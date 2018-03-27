#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=GenomeAnalysisTK
VERSION=4.0.2.1
#ARCHIVE=${SOFTWARE}-${VERSION}.tar.bz2
ARCHIVE=gatk-${VERSION}.zip
#echo "Prior to install the gatk module, you must download the archive $ARCHIVE manually, if not done already, from http://www.broadinstitute.org/gatk/download since it requires a license agreement.
#Once downloaded, copy it in \$MUGQIC_INSTALL_HOME_DEV/archive/ or \$MUGQIC_INSTALL_HOME/archive/"
SOFTWARE_DIR=${SOFTWARE}-${VERSION}
#ARCHIVE_URL=https://www.broadinstitute.org/gatk/download/auth?package=GATK
ARCHIVE_URL=https://github.com/broadinstitute/gatk/releases/download/${VERSION}/gatk-${VERSION}.zip
#JAR=$SOFTWARE.jar
JAR=gatk-package-${VERSION}-local.jar


# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
#  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
#  tar -jxvf $INSTALL_DOWNLOAD/$ARCHIVE --directory=$INSTALL_DIR/$SOFTWARE_DIR
  mkdir -p $INSTALL_DIR
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE -d $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend_path    PATH	            \$root
setenv          GATK_JAR            \$root/$JAR
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
