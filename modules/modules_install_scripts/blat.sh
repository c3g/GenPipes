

###################
################### BLAT
###################
VERSION="35x1"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/blat/blat-$VERSION
mkdir -p $INSTALL_PATH
wget "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat"
chmod +x blat
mv blat $INSTALL_PATH

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




