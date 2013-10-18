
###################
################### RSEM (http://deweylab.biostat.wisc.edu/rsem/)
###################
VERSION="1.2.6"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/rsem/
mkdir -p $INSTALL_PATH
cd  $INSTALL_PATH

# Download and extract
wget http://deweylab.biostat.wisc.edu/rsem/src/rsem-$VERSION.tar.gz
tar -xvf rsem-$VERSION.tar.gz
rm rsem-$VERSION.tar.gz 

# Compile
cd rsem-$VERSION
make -j8

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - RSEM (RNA-Seq by Expectation-Maximization) \"
}
module-whatis \"MUGQIC - RSEM \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/rsem/rsem-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/rsem
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/rsem/



