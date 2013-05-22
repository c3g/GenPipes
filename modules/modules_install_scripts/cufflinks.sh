
###################
################### Cufflinks
###################
VERSION="2.0.2"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/cufflinks
mkdir -p $INSTALL_PATH
wget http://cufflinks.cbcb.umd.edu/downloads/cufflinks-$VERSION.Linux_x86_64.tar.gz
tar -xvf cufflinks-$VERSION.Linux_x86_64.tar.gz
mv cufflinks-$VERSION.Linux_x86_64 $INSTALL_PATH


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Cufflinks \"
}
module-whatis \"MUGQIC - Cufflinks \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/cufflinks/cufflinks-$VERSION.Linux_x86_64
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/cufflinks
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/cufflinks/




