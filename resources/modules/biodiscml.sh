#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=BioDiscML
VERSION=1.9
ARCHIVE=${SOFTWARE}-${VERSION}.zip
ARCHIVE_URL=https://github.com/mickaelleclercq/BioDiscML/releases/download/v${VERSION}/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

build() {
  cd $INSTALL_DOWNLOAD

  # The release zip extracts flat (biodiscml.jar, config.conf, LICENSE.txt,
  # readme.md at its root, no wrapping folder), so unzip into a fresh
  # SOFTWARE_DIR ourselves before installing it.
  mkdir -p $SOFTWARE_DIR
  unzip $ARCHIVE -d $SOFTWARE_DIR

  # Install software
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE : machine-learning-based biomarker discovery (https://github.com/mickaelleclercq/BioDiscML)\"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          BIODISCML_HOME      \$root
setenv          BIODISCML_JAR       \$root/biodiscml.jar
prereq          mugqic/java/openjdk-jdk1.8.0_72
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
