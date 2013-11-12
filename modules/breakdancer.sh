

###################
################### Samtools
###################
VERSION="1.1.2_2013_03_08"

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/breakdancer/
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/breakdancer/tmp/unzip
cd $MUGQIC_INSTALL_HOME/modulefiles/mugqic/breakdancer/tmp

# Download

wget http://downloads.sourceforge.net/project/breakdancer/breakdancer-${VERSION}.zip
unzip breakdancer-${VERSION}.zip -d unzip/
mv unzip/breakdancer* unzip/breakdancer-${VERSION}
cd unzip/breakdancer-${VERSION}/

wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.6/samtools-0.1.6.tar.bz2
tar xvjf samtools-0.1.6.tar.bz2

# Compile samtools
cd samtools-0.1.6
make -j8

#install
cd ../breakdancer-${VERSION}/cpp
mv Makefile originalMakefile
sed "s|/Users/kchen3/pkg/samtools/samtools-0\.1\.6|\.\./samtools-0\.1\.6|g" originalMakefile > Makefile

make

cd ../../..

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/breakdancer # where to install..
ARCHIVE_PATH=$MUGQIC_INSTALL_HOME/archive/breakdancer 
mkdir -p $INSTALL_PATH $ARCHIVE_PATH
cp -r unzip/breakdancer-${VERSION}  $INSTALL_PATH
chmod -R 775 $INSTALL_PATH 
mv breakdancer-${VERSION}.zip $ARCHIVE_PATH
mv unzip/breakdancer-${VERSION}/samtools-0.1.6.tar.bz2 $ARCHIVE_PATH


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - reakdancer Structural Variant analyser \"
}
module-whatis \"Breakdancer Structural Variant analyser\"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/breakdancer/breakdancer-${VERSION}
prepend-path    PATH               \$root/cpp/
prepend-path    PATH               \$root/perl/
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version


mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/breakdancer/

cd ..
rm -rf tmp