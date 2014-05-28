#!/bin/bash

#
# BLAT
#

SOFTWARE=blat
VERSION=34

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV

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
# Write here the specific commands to download, extract, build the software, typically similar to:
ARCHIVE=${SOFTWARE}Src$VERSION.zip
wget http://users.soe.ucsc.edu/~kent/src/$ARCHIVE
unzip $ARCHIVE

SOFTWARE_DIR=$SOFTWARE-$VERSION
mkdir $INSTALL_DIR/$SOFTWARE_DIR
cd ${SOFTWARE}Src
# Won't compile elsewise
sed -i 's/ -Werror//g' inc/common.mk
export MACHTYPE=x86_64
make BINDIR=$INSTALL_DIR/$SOFTWARE_DIR

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX . $INSTALL_DIR/$SOFTWARE_DIR
mv -i $ARCHIVE ${!INSTALL_HOME}/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root
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
