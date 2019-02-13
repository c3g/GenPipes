t immediately on error
set -eu -o pipefail

SOFTWARE=bowtie
VERSION=1.2.0
#ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE=$SOFTWARE-$VERSION.zip
#ARCHIVE_URL=https://github.com/BenLangmead/${SOFTWARE}/archive/v${VERSION}.tar.gz
#ARCHIVE_URL=https://github.com/BenLangmead/${SOFTWARE}/releases/download/v${VERSION}.0/${SOFTWARE}-${VERSION}-src.zip
ARCHIVE_URL=https://sourceforge.net/projects/${SOFTWARE}-bio/files/${SOFTWARE}/${VERSION}/${SOFTWARE}-${VERSION%??}-source.zip
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
#  tar xzvf $ARCHIVE
  unzip $ARCHIVE

  mv ${SOFTWARE}-${VERSION%??} $SOFTWARE_DIR
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
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
