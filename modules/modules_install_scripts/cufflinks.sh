
###################
################### Cufflinks
###################
VERSION="2.1.1"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/cufflinks
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download and extract
wget http://cufflinks.cbcb.umd.edu/downloads/cufflinks-$VERSION.Linux_x86_64.tar.gz
tar zxvf cufflinks-$VERSION.Linux_x86_64.tar.gz
chmod -R g+w cufflinks-$VERSION.Linux_x86_64

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
rm cufflinks-$VERSION.Linux_x86_64.tar.gz
