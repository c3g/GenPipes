#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=java

# JDK 1.8
VERSION=openjdk-jdk1.8.0_72
ARCHIVE=jdk-8u72-ea-bin-b05-linux-x64-26_oct_2015.tar.gz
ARCHIVE_URL=http://download.java.net/jdk8u72/archive/b05/binaries/$ARCHIVE

# JDK 1.7
#VERSION=openjdk-jdk1.7.0_60
#ARCHIVE=jdk-7u60-ea-bin-b07-linux-x64-19_feb_2014.tar.gz
#ARCHIVE_URL=http://download.java.net/jdk7u60/archive/b07/binaries/$ARCHIVE

# JDK 1.6
#VERSION=openjdk-jdk1.6.0_38
#ARCHIVE=jdk-6u38-ea-bin-b04-linux-amd64-31_oct_2012.bin
#ARCHIVE_URL=http://download.java.net/jdk6/6u38/promoted/b04/binaries/$ARCHIVE

SOFTWARE_DIR=$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  # JDK 1.7 and 1.8
  tar zxvf $ARCHIVE
  # JDK 1.6
  #sh $ARCHIVE

  # Install software
  mv -i ${SOFTWARE_DIR/openjdk-/} $INSTALL_DIR/$SOFTWARE_DIR
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
