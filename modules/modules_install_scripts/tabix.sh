

###################
################### samtools
###################
VERSION="0.2.6"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/tabix/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download
wget http://downloads.sourceforge.net/project/samtools/tabix/tabix-$VERSION.tar.bz2
tar xvjf tabix-$VERSION.tar.bz2
# Compile
cd tabix-$VERSION
make -j8
# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - tabix \"
}
module-whatis \"Tabix indexer   \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/tabix/tabix-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tabix
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tabix/


