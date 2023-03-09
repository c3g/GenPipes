#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=purple
VERSION=3.7.1
ARCHIVE=${SOFTWARE}-${VERSION}.jar
ARCHIVE_URL=https://github.com/hartwigmedical/hmftools/releases/download/${SOFTWARE}-v${VERSION%.*}/${SOFTWARE}_v${VERSION}.jar
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD

  # Install software
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  cp $ARCHIVE $INSTALL_DIR/$SOFTWARE_DIR/${SOFTWARE}.jar
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE - require JDK1.8\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          PURPLE_HOME         \$root
setenv          PURPLE_JAR          \$root/purple.jar
setenv          PURPLE_DATA         /cvmfs/ref.mugqic/genomes/hmftools
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

