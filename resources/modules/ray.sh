#!/bin/bash

#
# Ray
#

SOFTWARE=ray
VERSION=2.3.1
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget http://downloads.sourceforge.net/project/denovoassembler/Ray-${VERSION}.tar.bz2
tar xjvf Ray-$VERSION.tar.bz2
cd Ray-$VERSION

# Guillimin
# module load gcc/4.7.2 openmpi/1.6.3-gcc
# Mammouth
# module load gcc/4.7.0 openmpi_gcc64/1.6.4
# Abacus
# module load itgenome/openmpi/1.6

make -j12 HAVE_LIBZ=y HAVE_LIBBZ2=y MAXKMERLENGTH=96 PREFIX=${INSTALL_PATH}/$SOFTWARE-$VERSION
make install
chmod -R a+rX,g+w ${INSTALL_PATH}/$SOFTWARE-$VERSION

# Add permissions and install software
cd $INSTALL_DOWNLOAD
mv -i Ray-$VERSION.tar.bz2 $MUGQIC_INSTALL_HOME_DEV/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;

set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$SOFTWARE-$VERSION
prepend-path    PATH                \$root;
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



