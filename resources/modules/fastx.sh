#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

###################
################### fastx
###################
# /sb/programs/analyste/software/fastx_toolkit-0.0.13.2
# http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit-0.0.13.2.tar.bz2
# http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2
# http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13.2_binaries_Linux_2.6_amd64.tar.bz2
VERSION="0.0.14"
LIBVERSION="0.7"
NAME=fastx_toolkit-$VERSION
LIBNAME=libgtextutils-$LIBVERSION
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/fastx/$NAME
INSTALL_DOWNLOAD=$MUGQIC_INSTALL_HOME/software/fastx/tmp
mkdir -p $INSTALL_PATH $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

wget https://github.com/agordon/fastx_toolkit/releases/download/$VERSION/$NAME.tar.bz2
wget https://github.com/agordon/libgtextutils/releases/download/$LIBVERSION/$LIBNAME.tar.gz
tar -xvf $NAME.tar.bz2
tar -xvf $LIBNAME.tar.gz

cd $LIBNAME
./configure --prefix=$INSTALL_PATH
make -j12
make install
cd ..

export PKG_CONFIG_PATH=$INSTALL_PATH/lib/pkgconfig
cd $NAME
./configure --prefix=$INSTALL_PATH
make -j12
make install
cd ..

chmod -R g+w $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - fastx \"
}
module-whatis \"MUGQIC - fastx \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/fastx/fastx_toolkit-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/fastx
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/fastx/
rm -rf $INSTALL_DOWNLOAD

