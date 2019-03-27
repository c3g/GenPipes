
###################
################### StringTie
###################
VERSION="1.3.4c"
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/stringtie
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download and extract
wget http://stringtie.cbcb.umd.edu/downloads/stringtie-$VERSION.Linux_x86_64.tar.gz
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-$VERSION.Linux_x86_64.tar.gz
tar zxvf stringtie-$VERSION.Linux_x86_64.tar.gz
chmod -R g+w stringtie-$VERSION.Linux_x86_64

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - StringTie \"
}
module-whatis \"MUGQIC - StringTie \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME_DEV)/software/stringtie/stringtie-$VERSION.Linux_x86_64
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/stringtie
mv .version $VERSION $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/stringtie/
rm stringtie-$VERSION.Linux_x86_64.tar.gz
