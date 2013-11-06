#!/bin/sh

#
# BVATools
#

SOFTWARE=bvatools
VERSION=1.1
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
cd $MUGQIC_INSTALL_HOME/archive

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget https://bitbucket.org/mugqic/${SOFTWARE}/downloads/${SOFTWARE}-${VERSION}.zip

mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
unzip $MUGQIC_INSTALL_HOME/archive/$SOFTWARE-$VERSION.zip

# Add permissions and install software
chmod -R ug+w .
chmod -R a+rX .

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE BAM and Variant analysis tools\" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;
setenv          BVATOOLS_HOME    \$root ;
setenv          BVATOOLS_JAR     \$root/$SOFTWARE-$VERSION-full.jar ;
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
