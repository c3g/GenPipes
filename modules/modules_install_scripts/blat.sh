

###################
################### BLAT
###################
VERSION="35"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/blat/blat-$VERSION/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
wget http://users.soe.ucsc.edu/~kent/src/blatSrc${VERSION}.zip
unzip blatSrc${VERSION}.zip
cd blatSrc
# Won't compile elsewise
sed -i 's/ -Werror//g' inc/common.mk
export MACHTYPE=x86_64
make BINDIR=${INSTALL_PATH}
cd ..
rm -r blatSrc

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - BLAT \"
}
module-whatis \"MUGQIC - BLAT \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/blat/blat-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/blat
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/blat/
