

###################
################### Trimmomatic
###################
VERSION="0.30"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/trimmomatic
mkdir -p $INSTALL_PATH
cd ${INSTALL_PATH}
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${VERSION}.zip
unzip Trimmomatic-$VERSION.zip

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - BEDtoolsTrimmomatic to trim fastqs \"
}
module-whatis \"Trimmomatic to trim fastqs  \"
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/trimmomatic/Trimmomatic-$VERSION
setenv          TRIMMOMATIC_JAR     \$root/trimmomatic-$VERSION.jar

" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trimmomatic
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trimmomatic/



