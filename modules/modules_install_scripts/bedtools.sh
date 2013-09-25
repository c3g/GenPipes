
###################
################### BEDTools
###################
VERSION="2.17.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/bedtools
INSTALL_DOWNLOAD=$MUGQIC_INSTALL_HOME/software/bedtools/tmp

mkdir -p $INSTALL_PATH $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download and extract
wget "http://bedtools.googlecode.com/files/BEDTools.v$VERSION.tar.gz"
tar -xvf BEDTools.v$VERSION.tar.gz

# Change the program directory name since different tar.gz archive versions have different names
mv *-$VERSION $INSTALL_PATH/bedtools-$VERSION
cd $INSTALL_PATH/bedtools-$VERSION
make -j8
cd ..
chmod -R g+w $INSTALL_PATH/bedtools-$VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - BEDtools \"
}
module-whatis \"MUGQIC - BEDtools  \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/bedtools/bedtools-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bedtools
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bedtools/
rm -rf $INSTALL_DOWNLOAD
