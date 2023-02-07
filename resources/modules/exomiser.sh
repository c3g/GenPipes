#!/bin/bash
# Exit immediately on error
set -eu -o pipefail


SOFTWARE=Exomiser
VERSION=7.2.1
ARCHIVE=${SOFTWARE,}-cli-${VERSION}-distribution.zip
ARCHIVE_URL=https://data.monarchinitiative.org/${SOFTWARE,}/${SOFTWARE,}-cli-${VERSION}-distribution.zip
SOFTWARE_DIR=${SOFTWARE}-${VERSION}

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD

build() {   
  cd $INSTALL_DOWNLOAD
  # Download the data (this is ~20GB and will take a while)
  wget -c https://data.monarchinitiative.org/${SOFTWARE,}/${SOFTWARE,}-cli-${VERSION}-data.zip

  # Download the checksums and verify the files (optional)
  wget -c https://data.monarchinitiative.org/${SOFTWARE,}/${SOFTWARE,}-cli-${VERSION}.sha256    
  sha256sum -c ${SOFTWARE,}-cli-${VERSION}.sha256

  # Unzip the distribution and data files - this will create a directory called 'exomiser-cli-7.2.1' in the current working directory
  unzip $ARCHIVE
  unzip ${SOFTWARE,}-cli-${VERSION}-data.zip
  mv  ${SOFTWARE,}-cli-$VERSION $SOFTWARE_DIR

  # Run a test genomiser analysis
  module_java=mugqic/java/openjdk-jdk1.8.0_72
  module load $module_java
  cd $SOFTWARE_DIR
  ln -s ${SOFTWARE,}-cli-${VERSION}.jar ${SOFTWARE,}.jar
  java -Xms2g -Xmx4g -jar ${SOFTWARE,}.jar --analysis NA19722_601952_AUTOSOMAL_RECESSIVE_POMP_13_29233225_5UTR_38.yml

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

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
setenv		EXOMISER_JAR        \$root/${SOFTWARE,}.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
