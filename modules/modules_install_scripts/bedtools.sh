
###################
################### BEDTools
###################
VERSION="2.17.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/bedtools/bedtools-$VERSION
mkdir -p $INSTALL_PATH
wget "http://bedtools.googlecode.com/files/BEDTools.v$VERSION.tar.gz"
tar -xvf BEDTools.v$VERSION.tar.gz
cd bedtools-$VERSION
make -j8
mv ./* $INSTALL_PATH



# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - BEDtools \"
}
module-whatis \"MUGQIC - BEDtools  \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/bedtools/bedtools-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bedtools
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bedtools/





