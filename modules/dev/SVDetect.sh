

######### VirusFinder ##################
VERSION="0.8b"
NAME="SVDetect" # same could apply to all ucsc tools
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/sv_detect/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download and extract 
wget http://sourceforge.net/projects/svdetect/files/latest/download?source=files
tar -xvzf ${NAME}_r$VERSION.tar.gz

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - SVDetect \"
}
module-whatis \"MUGQIC - SVDetect \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/sv_detect/SVDetect_r$VERSION/bin/
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/svdetect
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/svdetect/
