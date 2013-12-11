#!/bin/sh

#
# RSEM
#

SOFTWARE=rsem
VERSION=1.2.8
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget http://deweylab.biostat.wisc.edu/rsem/src/$SOFTWARE-$VERSION.tar.gz
tar zxvf $SOFTWARE-$VERSION.tar.gz
cd $SOFTWARE-$VERSION
make

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
chmod -R o+rX .
mv -i $SOFTWARE-$VERSION $INSTALL_PATH
mv -i $SOFTWARE-$VERSION.tar.gz $MUGQIC_INSTALL_HOME_DEV/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$SOFTWARE-$VERSION ;
prepend-path    PATH                \$root ;
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
