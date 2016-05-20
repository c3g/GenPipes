#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


# "http://topaz.gatech.edu/GeneMark/tmp/GMtool_XMhCb/gm_et_linux_64.tar.gz"
# "http://topaz.gatech.edu/GeneMark/tmp/GMtool_XMhCb/gm_key_64.gz"
#
#
# "http://topaz.gatech.edu/GeneMark/tmp/GMtool_y9U2C/genemark_suite_linux_64.tar.gz"
# "http://topaz.gatech.edu/GeneMark/tmp/GMtool_y9U2C/gm_key_64.gz"

# rm -rf /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_dev/software/genemark-es/genemark-es-4.32  /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_dev/modulefiles/mugqic_dev/genemark-es/4.32


################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE="genemark-es"
VERSION="4.32"
ARCHIVE="gm_et_linux_64.tar.gz"
ARCHIVE_URL="http://topaz.gatech.edu/GeneMark/tmp/GMtool_XMhCb/$ARCHIVE"
SOFTWARE_DIR=$SOFTWARE-$VERSION 

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC COMMAND
  mv gm_et_linux_64 $SOFTWARE_DIR

  wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_XMhCb/gm_key_64.gz -O gm_key_64.gz
	gunzip gm_key_64.gz
	mv gm_key_64  $SOFTWARE_DIR/

  # Install software
  cd $INSTALL_DOWNLOAD  ## TO BE ADDED AND MODIFIED IF NECESSARY
  mv -i $SOFTWARE_DIR $INSTALL_DIR/  ## TO BE ADDED AND MODIFIED IF NECESSARY
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"
set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/gmes_petap
set            home                 \$::env(HOME)
file delete \$home/.gm_key
file copy       $INSTALL_DIR/$SOFTWARE_DIR/gmes_petap/gm_key \$home/.gm_key
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
