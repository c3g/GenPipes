#!/bin/bash

#
# wgs-assembler
#
SOFTWARE=wgs-assembler
VERSION=8.1
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE/
# Don't use /tmp!!! we need to compile wgs-assembler in it's final location because it hard codes the path in sources
mkdir -p $INSTALL_PATH
cd $MUGQIC_INSTALL_HOME/archive/

# Download, extract, build
wget http://downloads.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-${VERSION}/wgs-${VERSION}.tar.bz2
cd $INSTALL_PATH
tar xjvf $MUGQIC_INSTALL_HOME/archive/wgs-${VERSION}.tar.bz2
cd wgs-${VERSION}/kmer
make install
cd ../samtools
make
cd ../src
# change default size to 131k. Default since 8.1 is 65k
sed -i 's/\(#define AS_READ_MAX_NORMAL_LEN_BITS[^0-9]\+\)[0-9]\+/\1 17/g' AS_global.H
make

# Add permissions and install software
cd $INSTALL_PATH
chmod -R a+rX .
chmod -R ug+w .

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/wgs-${VERSION}
prepend-path    PATH                \$root/Linux-amd64/bin/
" > $VERSION

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
