#TODO:
# Support Abacus, Mammouth


###################
################### ABYSS
################### 
VERSION="1.3.5"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/abyss/abyss-$VERSION
mkdir -p $INSTALL_PATH
wget "http://www.bcgsc.ca/platform/bioinfo/software/abyss/releases/$VERSION/abyss-$VERSION.tar.gz"
mv abyss-$VERSION.tar.gz $INSTALL_PATH
cd $INSTALL_PATH
tar -xvzf $INSTALL_PATH/abyss-$VERSION.tar.gz


#guillimin
module add gcc/4.7.2
module load openmpi/1.6.3-gcc
PATH_OPENMPI_LIB="/software/CentOS-5/tools/openmpi-1.6.3-gcc"
#mammouth : TODO error with mpi in configure step
#module add gcc/4.7.0
#module load openmpi_gcc64/1.6.4
PATH_OPENMPI_LIB="/opt/mpi/gcc/openmpi-1.6.4"
#abacus : TODO error with mpi in configure step 
#module add compat-openmpi-x86_64
#PATH_OPENMPI_LIB="/usr/lib64/compat-openmpi/lib"

#TODO : for other cluster

cd $INSTALL_PATH/abyss-$VERSION

#Maybe you need to change this part for new version
wget http://downloads.sourceforge.net/project/boost/boost/1.50.0/boost_1_50_0.tar.bz2
tar jxf boost_1_50_0.tar.bz2
ln -s boost_1_50_0/boost boost

./configure --prefix=$INSTALL_PATH --with-mpi=$PATH_OPENMPI_LIB
make
make install

cd ~

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - ABySS \"
}
module-whatis \"MUGQIC - ABySS \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/abyss/abyss-$VERSION/bin/
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/abyss
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/abyss/

###################
################### ray 
###################

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/ray/
VERSION="2.2.0-rc0"

cd $INSTALL_PATH/ray-$VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - RAY \"
}
module-whatis \"MUGQIC - RAY \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/ray/ray-$VERSION/
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/ray
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/ray/

