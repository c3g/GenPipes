#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=picard
#version 2 or later require JDK1.8
VERSION=2.17.3
#ARCHIVE=$SOFTWARE-tools-$VERSION.zip								# for version < 2.6.0
#ARCHIVE_URL=https://github.com/broadinstitute/picard/releases/download/$VERSION/$ARCHIVE       # for version < 2.6.0
#SOFTWARE_DIR=$SOFTWARE-tools-$VERSION								# for version < 2.6.0
ARCHIVE=${SOFTWARE}-${VERSION}.jar								# for version > 2.5.0
ARCHIVE_URL=https://github.com/broadinstitute/picard/releases/download/$VERSION/${SOFTWARE}.jar	# for version > 2.5.0
SOFTWARE_DIR=$SOFTWARE-$VERSION									# for version > 2.5.0

build() {
  cd $INSTALL_DOWNLOAD
#  unzip $ARCHIVE										# for version < 2.6.0

  # Install software
#  mv -i $SOFTWARE_DIR $INSTALL_DIR/                                                            # for version < 2.6.0
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR								# for version > 2.5.0
  cp $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/${SOFTWARE}.jar					# for version > 2.5.0
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE - require JDK1.8\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          PICARD_HOME         \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

