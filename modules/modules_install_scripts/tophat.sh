
###################
################### Tophat
###################
VERSION="2.0.9"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/tophat
INSTALL_DOWNLOAD=$MUGQIC_INSTALL_HOME/software/tophat/tmp
mkdir -p $INSTALL_PATH $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download and extract
wget http://tophat.cbcb.umd.edu/downloads/tophat-$VERSION.Linux_x86_64.tar.gz
tar -xvf tophat-$VERSION.Linux_x86_64.tar.gz

# Install
mv tophat-$VERSION.Linux_x86_64/ $INSTALL_PATH
cd $INSTALL_PATH
chmod -R g+w tophat-$VERSION.Linux_x86_64/

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Adds tophat to your environment \"
}
module-whatis \"MUGQIC - Adds tophat to your environment \"
                       
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/tophat/tophat-$VERSION.Linux_x86_64
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tophat
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tophat
rm -rf $INSTALL_DOWNLOAD
