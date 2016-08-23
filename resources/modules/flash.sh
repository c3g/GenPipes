#!/bin/sh
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=FLASH
VERSION=1.2.8
INSTALL_PATH=$MUGQIC_INSTALL_HOME_TMP/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget https://sourceforge.net/projects/flashpage/files/$SOFTWARE-$VERSION.tar.gz
tar -zxvf $SOFTWARE-$VERSION.tar.gz                                                  
cd $SOFTWARE-$VERSION                                                            
make                                                                             

# Add permissions and install software
chmod -R 775 *
cd $INSTALL_DOWNLOAD
mv -i $SOFTWARE-$VERSION $INSTALL_PATH                                          
mv -i $INSTALL_DOWNLOAD/$SOFTWARE-$VERSION.tar.gz $MUGQIC_INSTALL_HOME_TMP/archive      

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE-$VERSION \" ;
}
module-whatis \"$SOFTWARE  \" ; 
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_TMP)/software/$SOFTWARE/$SOFTWARE-$VERSION ;
prepend-path    PATH                \$root/bin ;  
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_TMP/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_TMP/modulefiles/mugqic/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
