#!/bin/bash

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


#
# Software_name repeatmasker
#

SOFTWARE=RepeatMasker
VERSION=open-4-0-5
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget http://www.repeatmasker.org/RepeatMasker-$VERSION.tar.gz
tar -xvf $SOFTWARE-$VERSION.tar.gz                                             
cd $SOFTWARE                                                            
./configure --prefix=$INSTALL_PATH/$SOFTWARE-$VERSION                            
make                                                                             
make install

# Add permissions and install software
chmod -R 775 *
cd $INSTALL_DOWNLOAD
#mv -i $SOFTWARE-$VERSION $INSTALL_PATH                                          
mv -i $INSTALL_DOWNLOAD/$SOFTWARE-$VERSION.tar.gz $MUGQIC_INSTALL_HOME_DEV/archive      

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE-$VERSION \" ;
}
module-whatis \"$SOFTWARE  \" ; 
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$SOFTWARE-$VERSION ;
prepend-path    PATH                \$root/bin ;  
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
