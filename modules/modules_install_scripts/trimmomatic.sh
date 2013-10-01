

###################
################### Trimmomatic
###################
VERSION=0.30
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/trimmomatic
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_PATH/archive $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download and extract
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$VERSION.zip
unzip Trimmomatic-$VERSION.zip

# Install
mv Trimmomatic-$VERSION $INSTALL_PATH
cd $INSTALL_PATH
chmod -R g+w Trimmomatic-$VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Trimmomatic to trim fastq \"
}
module-whatis \"Trimmomatic to trim fastq  \"
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/trimmomatic/Trimmomatic-$VERSION
setenv          TRIMMOMATIC_JAR     \$root/trimmomatic-$VERSION.jar

" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trimmomatic
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trimmomatic/
mv $INSTALL_DOWNLOAD/Trimmomatic-$VERSION.zip $INSTALL_PATH/archive/
rm -rf $INSTALL_DOWNLOAD
