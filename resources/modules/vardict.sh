#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

HOST=`hostname`;
DNSDOMAIN=`dnsdomainname`;

SOFTWARE=VarDictJava
VERSION=1.8.2
ARCHIVE=$SOFTWARE-$VERSION.zip
ARCHIVE_URL=https://github.com/AstraZeneca-NGS/${SOFTWARE}/releases/download/v${VERSION}/VarDict-${VERSION}.zip
SOFTWARE_DIR=$SOFTWARE-$VERSION 

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
  mv -i VarDict-${VERSION}  $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin;
setenv          VARDICT_HOME        \$root;
setenv          VARDICT_BIN         \$root/bin;
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
