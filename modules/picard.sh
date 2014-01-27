#!/bin/sh

#
# Picard
#

SOFTWARE=picard
VERSION=1.106
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget http://downloads.sourceforge.net/project/picard/picard-tools/$VERSION/$SOFTWARE-tools-$VERSION.zip
unzip $SOFTWARE-tools-$VERSION.zip

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
mv -i $SOFTWARE-tools-$VERSION $INSTALL_PATH
mv -i $SOFTWARE-tools-$VERSION.zip $MUGQIC_INSTALL_HOME/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-tools-$VERSION ;
setenv          PICARD_HOME         \$root ;
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
