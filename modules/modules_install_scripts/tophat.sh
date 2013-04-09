
###################
################### tophat
###################
VERSION="2.0.8"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/tophat/tophat-$VERSION.Linux_x86_64 # where to install..
mkdir -p $INSTALL_PATH
wget http://tophat.cbcb.umd.edu/downloads/tophat-$VERSION.Linux_x86_64.tar.gz
tar -xvf tophat-$VERSION.Linux_x86_64.tar.gz
mv tophat-$VERSION.Linux_x86_64/* $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Adds tophat to your environment \"
}
module-whatis \"MUGQIC - Adds tophat to your environment \"
                       
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/tophat/tophat-$VERSION.Linux_x86_64
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tophat
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tophat


