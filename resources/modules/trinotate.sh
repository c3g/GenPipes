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
VERSION=20131110
ARCHIVE=${SOFTWARE^}_r$VERSION.tar.gz
ARCHIVE_URL=http://downloads.sourceforge.net/project/$SOFTWARE/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE^}_r$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Download Trinotate Sqlite template DB
  SQLITE_ARCHIVE=Trinotate.sqlite.gz
  download_archive "http://sourceforge.net/projects/trinotate/files/TRINOTATE_RESOURCES/TrinotateSqlite.sprot.$VERSION.db.gz/download" $SQLITE_ARCHIVE
  gunzip $SQLITE_ARCHIVE -c > $SOFTWARE_DIR/${SQLITE_ARCHIVE/.gz/}

  # Install software
  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  store_archive $SQLITE_ARCHIVE
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
setenv          TRINOTATE_HOME      \$root
setenv          TRINOTATE_SQLITE    \$root/Trinotate.sqlite
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
