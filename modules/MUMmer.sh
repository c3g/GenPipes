#!/bin/sh

#
# Mummer
#
SOFTWARE=MUMmer
VERSION=3.23
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE/
# Don't use /tmp!!! we need to compile MUMmer in it's final location because it hard codes the path in sources
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
wget http://downloads.sourceforge.net/project/mummer/mummer/${VERSION}/${SOFTWARE}${VERSION}.tar.gz
tar xzvf ${SOFTWARE}${VERSION}.tar.gz
cd ${SOFTWARE}$VERSION
make check
make install

# Add permissions and install software
cd $INSTALL_PATH
chmod -R ug+rwX .
chmod -R o+rX .
mv -i ${SOFTWARE}${VERSION}.tar.gz $MUGQIC_INSTALL_HOME/archive/

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/${SOFTWARE}${VERSION}
prepend-path    PATH                \$root
prepend-path    PERL5LIB            \$root/scripts/
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
