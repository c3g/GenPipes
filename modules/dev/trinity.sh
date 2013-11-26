#!/bin/sh

#
# Trinity
#

SOFTWARE=trinity
VERSION=20131110
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget http://sourceforge.net/projects/trinityrnaseq/files/trinityrnaseq_r$VERSION.tar.gz
tar zxvf trinityrnaseq_r$VERSION.tar.gz
cd trinityrnaseq_r$VERSION
make

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
mv -i trinityrnaseq_r$VERSION $INSTALL_PATH
mv -i trinityrnaseq_r$VERSION.tar.gz $MUGQIC_INSTALL_HOME_DEV/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/trinityrnaseq_r$VERSION ;
setenv          TRINITY_HOME        \$root ;
prepend-path    PATH                \$root ;
prepend-path    PATH                \$root/util ;
prepend-path    PATH                \$root/util/RSEM_util ;
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
chmod -R ug+rwX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
