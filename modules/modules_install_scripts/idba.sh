


###################
################### IDBA
###################
VERSION="1.1.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/idba/
mkdir -p $INSTALL_PATH


# IDBA : compile and test whether it segfaults
wget http://hku-idba.googlecode.com/files/idba-$VERSION.tar.gz
tar -xvf idba-$VERSION.tar.gz
cd idba-$VERSION



./configure --prefix=$INSTALL_PATH/idba-$VERSION
make -j8
make install
cp -rf bin/* $INSTALL_PATH/idba-$VERSION/bin


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - IDBA assembler \"
}
module-whatis \"MUGQIC - IDBA assembler  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/idba/idba-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/idba
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/idba/






###################
################### IDBA LONG
###################
VERSION="1.0.9"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/idba_long/
mkdir -p $INSTALL_PATH


# IDBA : compile and test whether it segfaults
wget http://hku-idba.googlecode.com/files/idba_ud-$VERSION.tar.gz
tar -xvf idba_ud-$VERSION.tar.gz
cd idba_ud-$VERSION

#wget http://hku-idba.googlecode.com/files/idba-$VERSION.tar.gz
#tar -xvf idba-$VERSION.tar.gz
#cd idba-$VERSION

# TODO MANUAL: Tweak for longer reads https://groups.google.com/forum/?fromgroups=#!topic/hku-idba/ShB95FPswN8
# https://groups.google.com/forum/?fromgroups=#!topic/hku-idba/NE2JXqNvTFY


./configure --prefix=$INSTALL_PATH/idba-$VERSION
make -j8
make install
cp -rf bin/* $INSTALL_PATH/idba-$VERSION/bin


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - IDBA assembler \"
}
module-whatis \"MUGQIC - IDBA assembler  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/idba_long/idba-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/idba_long
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/idba_long/




