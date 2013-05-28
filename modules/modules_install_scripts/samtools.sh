

###################
################### samtools
###################
VERSION="0.1.19"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/samtools/
mkdir -p $INSTALL_PATH
# Download
wget http://downloads.sourceforge.net/project/samtools/samtools/$VERSION/samtools-$VERSION.tar.bz2
tar xvjf samtools-$VERSION.tar.bz2
# Compile
cd samtools-$VERSION
make -j8
cd ..
# move
mv samtools-$VERSION $INSTALL_PATH
# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - samtools \"
}
module-whatis \"SAMTools SAM/BAM manipulation   \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/samtools/samtools-$VERSION
prepend-path    PATH               \$root
prepend-path    PATH               \$root/bcftools
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/samtools
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/samtools/


