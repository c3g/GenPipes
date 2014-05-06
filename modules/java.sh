#!/bin/bash

#
# Java
#

SOFTWARE=java
VERSION=openjdk-jdk1.7.0_60
#VERSION=openjdk-jdk1.6.0_38

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir $INSTALL_DOWNLOAD
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

SOFTWARE_DIR=$VERSION
# Prefix default "jdk..." software directory by "openjdk-" to differentiate it from oracle-jdk
mv ${SOFTWARE_DIR/openjdk-/} $SOFTWARE_DIR

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX .
mv -i $SOFTWARE_DIR $INSTALL_DIR
mv -i $ARCHIVE ${!INSTALL_HOME}/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE from OpenJDK \" ;
}
module-whatis \"$SOFTWARE from OpenJDK \" ;

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin ;
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by removing '_INSTALL_HOME' in $INSTALL_HOME and lowercasing the result
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME/_INSTALL_HOME/} | tr '[:upper:]' '[:lower:]'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX $VERSION .version
mv $VERSION .version $MODULE_DIR

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
