# Exit immediately on error
set -eu -o pipefail

SOFTWARE=bambam
VERSION=1.4
ARCHIVE=$SOFTWARE-${VERSION}.tar.gz
ARCHIVE_URL=https://sourceforge.net/projects/bambam/files/$SOFTWARE-${VERSION}.tgz
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar -xzvf $ARCHIVE

  mv $SOFTWARE $SOFTWARE_DIR
  cd $SOFTWARE_DIR
  make -j12 all

  cd $INSTALL_DOWNLOAD
  mv $SOFTWARE_DIR $INSTALL_DIR/

}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root        $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH        \$root/bin
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
