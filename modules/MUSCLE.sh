#!/bin/bash

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


#
# Software_name  gnuplot.
#

SOFTWARE=MUSCLE
VERSION=3.8.31
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE/$SOFTWARE-$VERSION
INSTALL_DOWNLOAD=$SCRATCH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget http://www.drive5.com/muscle/downloads${VERSION}/muscle${VERSION}_i86linux64.tar.gz
tar -xvf muscle${VERSION}_i86linux64.tar.gz   
mkdir -p $INSTALL_PATH/bin                                       
cp muscle${VERSION}_i86linux64 $INSTALL_PATH/bin/
chmod -R 775 $INSTALL_PATH

# Add permissions and install software
chmod -R 775 *

cd $INSTALL_PATH/bin
ln -s muscle${VERSION}_i86linux64 muscle 

cd $INSTALL_DOWNLOAD
mv -i $INSTALL_DOWNLOAD/muscle${VERSION}_i86linux64.tar.gz $MUGQIC_INSTALL_HOME/archive/ 

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE-$VERSION \" ;
}
module-whatis \"$SOFTWARE  \" ; 
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;
prepend-path    PATH                \$root/bin ;  
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
