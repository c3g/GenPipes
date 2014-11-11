#!/bin/bash

#
# Aspera
#

SOFTWARE=aspera-connect
VERSION=3.3.3
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE/$SOFTWARE-$VERSION/
mkdir -p ${INSTALL_PATH}

# Download, extract, build
echo "The installation needs to be done manually"
echo "Download aspera-connect http://downloads.asperasoft.com/connect2/"
echo "Then excute the script and it will installe in your $HOME/.aspera/connect"
echo "then rsync -avP $HOME/.aspera/connect/* $INSTALL_PATH"
echo ""
echo "This script already created the module"
echo "don't forget to chmod -R a+rX $INSTALL_PATH"

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;
prepend-path    PATH                \$root/bin ;  ## TO BE ADDED IF NECESSARY
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

