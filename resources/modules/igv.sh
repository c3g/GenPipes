#!/bin/bash

###################
################### IGV
###################
VERSION="2.3.23"
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/igv/tmp/unzip/
cd $MUGQIC_INSTALL_HOME/modulefiles/mugqic/igv/tmp

# Download and install
wget http://www.broadinstitute.org/igv/projects/downloads/IGV_${VERSION}.zip
unzip IGV_$VERSION.zip -d unzip/
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/igv # where to install..
ARCHIVE_PATH=$MUGQIC_INSTALL_HOME/archive/igv 
mkdir -p $INSTALL_PATH $ARCHIVE_PATH
cp -r unzip/IGV_${VERSION}  $INSTALL_PATH
chmod -R 775 $INSTALL_PATH 
mv IGV_${VERSION}.zip $ARCHIVE_PATH

# Module filem
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - IGV \"
}
module-whatis \"MUGQIC - IGVtools  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/igv/IGV_${VERSION}
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/igv/

cd ..
rm -rf tmp
