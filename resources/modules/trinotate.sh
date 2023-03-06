#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# NOTE: 
# - Perl module DB_File is a dependency for the Transdecoder part of Trinotate. This module depends on some version BerkeleyDB which was not present on Mammouth...
#
# NOTES:
# - Assuming trinotate and trinity version follow one another
# - Assuming sqlite is already available on the system

SOFTWARE=trinotate
VERSION=4.0.0
ARCHIVE=${SOFTWARE^}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/${SOFTWARE^}/${SOFTWARE^}/archive/refs/tags/${SOFTWARE^}-v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE^}-${VERSION}

MODULE_PERL=mugqic/perl/5.34.0
MODULE_INFERNAL=mugqic/infernal/1.1.4
MODULE_DIAMOND=mugqic/diamond/2.1.4


build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  cd ${SOFTWARE^}-${SOFTWARE^}-v${VERSION}
  module load $MODULE_PERL $MODULE_INFERNAL $MODULE_DIAMOND
  SQLITE=Trinotate.sqlite
  ./Trinotate \
    --create \
    --db ${SQLITE} \
    --trinotate_data_dir data \
    --use_diamond

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i ${SOFTWARE^}-${SOFTWARE^}-v${VERSION} $INSTALL_DIR/$SOFTWARE_DIR

  # Download Trinotate resources (adjust file names for newer Trinotate version)
  # SQLITE=Trinotate.sprot_uniref90.20150131.boilerplate.sqlite
  # download_archive "ftp://ftp.broadinstitute.org/pub/users/bhaas/Trinotate_v${VERSION%.*}_RESOURCES" $SQLITE.gz
  # gunzip $SQLITE.gz -c > $INSTALL_DIR/$SOFTWARE_DIR/$SQLITE
  # store_archive $SQLITE.gz
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - ${SOFTWARE^} \"
}
module-whatis \"${SOFTWARE^}\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
setenv          TRINOTATE_HOME      \$root
setenv          TRINOTATE_SQLITE    \$root/$SQLITE
setenv          TRINOTATE_DATA_DIR  \$root/data
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
