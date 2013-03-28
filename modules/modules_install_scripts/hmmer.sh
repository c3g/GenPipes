
###################
################### HMMER3
###################
VERSION="3.0"
INSTALL_PATH="$MUGQIC_INSTALL_HOME/software/hmmer/hmmer-"$VERSION
mkdir -p $INSTALL_PATH
FILE="hmmer-"$VERSION"-linux-intel-x86_64.tar.gz" 
wget "http://selab.janelia.org/software/hmmer3/$VERSION/$FILE" 
tar -xvf $FILE
cd "hmmer-$VERSION-linux-intel-x86_64"
./configure
make
make check
cd ..
mv "hmmer-$VERSION-linux-intel-x86_64" $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - hmmer3 \"
}
module-whatis \"MUGQIC - hmmer3 \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/hmmer/hmmer-$VERSION/hmmer-$VERSION-linux-intel-x86_64/binaries
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/hmmer
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/hmmer/




