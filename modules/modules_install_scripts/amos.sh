

###################
################### BWA
###################
# tpx patch can be found here:
# ftp://ftp.conveysupport.com/outgoing/bwa/bwa-0.6.2-tpx.patch
VERSION="3.1.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/amos/
mkdir -p $INSTALL_PATH

# Download
cd $MUGQIC_INSTALL_HOME/archive/
wget http://downloads.sourceforge.net/project/amos/amos/${VERSION}/amos-${VERSION}.tar.gz

cd $INSTALL_PATH
tar xvzf $MUGQIC_INSTALL_HOME/archive/amos-${VERSION}.tar.gz

# compile
mv amos-${VERSION} amos-${VERSION}-src
cd amos-${VERSION}-src
# with-qmake seems to work on guillimin and mammouth as well as abacus
module load mugqic/mummer mugqic/ucsc
./configure --prefix $MUGQIC_INSTALL_HOME/software/amos/amos-${VERSION} --with-qmake-qt4=/usr/lib64/qt4/bin/qmake
make
make install && \
chmod -R ug+rwX $MUGQIC_INSTALL_HOME/software/amos/amos-${VERSION} && \
rm -rf ../amos-${VERSION}-src

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Amos \"
}
module-whatis \"MUGQIC - Amos  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/amos/amos-${VERSION}
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/amos
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/amos/

