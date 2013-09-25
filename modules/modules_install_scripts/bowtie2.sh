###################
################### Bowtie2
###################
VERSION="2.1.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/bowtie/bowtie-$VERSION
INSTALL_DOWNLOAD=$MUGQIC_INSTALL_HOME/software/bowtie/tmp
mkdir -p $INSTALL_PATH $INSTALL_DOWNLOAD

# Download and extract
cd $INSTALL_DOWNLOAD
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$VERSION/bowtie2-$VERSION-source.zip/download
unzip bowtie2-$VERSION-source.zip

# Compile
cd bowtie2-$VERSION
make -j8
mv bowtie2* $INSTALL_PATH
chmod -R g+w $INSTALL_PATH
cd ..


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Bowtie2 aligner \"
}
module-whatis \"MUGQIC - Bowtie2 aligner \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/bowtie/bowtie-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bowtie
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bowtie/

rm -rf $INSTALL_DOWNLOAD
