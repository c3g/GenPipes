#!/bin/bash

#
# Java
#

SOFTWARE=java
VERSION=jdk1.7.0_60
#VERSION=jdk1.6.0_38
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build

# JDK 1.7
ARCHIVE=jdk-7u60-ea-bin-b07-linux-x64-19_feb_2014.tar.gz
wget http://download.java.net/jdk7u60/archive/b07/binaries/$ARCHIVE
tar zxvf $ARCHIVE

# JDK 1.6
#ARCHIVE=jdk-6u38-ea-bin-b04-linux-amd64-31_oct_2012.bin
#wget http://download.java.net/jdk6/6u38/promoted/b04/binaries/$ARCHIVE
#sh $ARCHIVE

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
chmod -R o+rX .
mv -i $VERSION $INSTALL_PATH
mv -i $ARCHIVE $MUGQIC_INSTALL_HOME/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE from OpenJDK \"
}
module-whatis \"$SOFTWARE from OpenJDK  \"

set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$VERSION
prepend-path    PATH                \$root/bin
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
