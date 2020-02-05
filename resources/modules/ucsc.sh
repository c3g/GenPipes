#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ucsc
# By default, the latest remote version will be downloaded and the version date set appropriately.
# To use a local archive specific version, uncomment and update VERSION
VERSION=v387
#VERSION=latest
#VERSION=20140212
if [[ $VERSION == latest ]]
then
  ARCHIVE=userApps.src.tgz
else
  if [[ ${VERSION:0:1} == "v" ]]
  then
    ARCHIVE=userApps.$VERSION.src.tgz
  else
    ARCHIVE=ucsc-userApps-$VERSION.src.tgz
  fi
fi
ARCHIVE_URL=http://hgdownload.cse.ucsc.edu/admin/exe/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD

  if [[ $VERSION == latest ]]
  then
    REMOTE_ARCHIVE=$ARCHIVE
    # Set VERSION with the archive last modification date
    VERSION=`stat --printf=%y $REMOTE_ARCHIVE | perl -pe 's/^(\d+)-(\d+)-(\d+).*/\1\2\3/'`
    ARCHIVE=$SOFTWARE-userApps-$VERSION.src.tgz
    mv $REMOTE_ARCHIVE $ARCHIVE
    SOFTWARE_DIR=$SOFTWARE-$VERSION
  fi

  tar zxvf $ARCHIVE
  mv userApps $SOFTWARE_DIR

  cd $SOFTWARE_DIR
  make -j12

  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
