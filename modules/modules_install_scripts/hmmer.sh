
###################
################### HMMER3
###################
VERSION=3.0
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/hmmer/hmmer-$VERSION
INSTALL_DOWNLOAD=$MUGQIC_INSTALL_HOME/software/hmmer/tmp
mkdir -p $INSTALL_PATH $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Adjust remote download URL according to version first number
if [[ ${VERSION:0:1} == 3 ]]
then
  SUFFIX=3
else
  SUFFIX=""
fi

# Download and extract
wget http://selab.janelia.org/software/hmmer$SUFFIX/$VERSION/hmmer-$VERSION.tar.gz 
tar -xvf hmmer-$VERSION.tar.gz

# Compile and install
cd hmmer-$VERSION
./configure --prefix $INSTALL_PATH
make
make check
make install
cd ..
chmod -R g+w $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - hmmer-$VERSION \"
}
module-whatis \"MUGQIC - hmmer-$VERSION \"

set             root               \$::env(MUGQIC_INSTALL_HOME)/software/hmmer/hmmer-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/hmmer
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/hmmer/
rm -rf $INSTALL_DOWNLOAD
