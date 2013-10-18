

###################
################### Ray
###################
#module load compat-openmpi-x86_64 # 4abacus
#module load compat-openmpi-psm-x86_64
VERSION="2.2.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/ray/ray-$VERSION-k96
wget http://downloads.sourceforge.net/project/denovoassembler/Ray-v$VERSION.tar.bz2
tar xjvf Ray-v$VERSION.tar.bz2
mv Ray-v$VERSION Ray-v$VERSION-src
cd Ray-v$VERSION-src

make HAVE_LIBZ=y HAVE_LIBBZ2=y PREFIX=${INSTALL_PATH} MAXKMERLENGTH=96
make HAVE_LIBZ=y HAVE_LIBBZ2=y PREFIX=${INSTALL_PATH} MAXKMERLENGTH=96 install

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Ray assembler \"
}
module-whatis \"Ray Parallel genome assemblies for parallel DNA sequencing \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/ray/ray-$VERSION-k96
prepend-path    PATH               \$root
" > $VERSION-k96

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION-k96\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/ray
mv .version $VERSION-k96 $MUGQIC_INSTALL_HOME/modulefiles/mugqic/ray/


