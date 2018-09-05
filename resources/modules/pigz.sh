#!/bin/bash

###################
################### exonerate
###################
VERSION="2.3.3"
INSTALL_PATH=$MUGQIC_INSTALL_HOME_TMP/software/pigz/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download
wget --content-disposition http://zlib.net/pigz/pigz-${VERSION}.tar.gz
tar xvzf pigz-${VERSION}.tar.gz
# Compile
cd pigz-${VERSION}
make -j12

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - pigz \"
}
module-whatis \"A parallel implementation of gzip for modern multi-processor, multi-core machines\"
            
set             root               \$::env(MUGQIC_INSTALL_HOME_TMP)/software/pigz/pigz-${VERSION}
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME_TMP/modulefiles/mugqic/pigz
mv .version $VERSION $MUGQIC_INSTALL_HOME_TMP/modulefiles/mugqic/pigz/


