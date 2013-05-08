

###################
################### exonerate
###################
VERSION="2.2.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/exonerate/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download
wget http://www.ebi.ac.uk/~guy/exonerate/exonerate-${VERSION}.tar.gz
tar xvzf exonerate-${VERSION}.tar.gz
# Compile
mv exonerate-${VERSION} exonerate-${VERSION}-src
cd exonerate-${VERSION}-src
./configure --prefix ${INSTALL_PATH}/exonerate-${VERSION}
make -j8
make install && ../exonerate-${VERSION}-src
# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - exonerate \"
}
module-whatis \"Exonerate is a generic tool for pairwise sequence comparison and fasta manipulation\"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/exonerate/exonerate-${VERSION}
prepend-path    PATH               \$root/bin
prepend-path    MANPATH            \$root/share/man/
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/exonerate
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/exonerate/


