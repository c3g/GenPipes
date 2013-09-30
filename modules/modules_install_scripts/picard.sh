
###################
################### Picard
###################
VERSION="1.98"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/picard/
INSTALL_DOWNLOAD=$MUGQIC_INSTALL_HOME/software/picard/tmp
mkdir -p $INSTALL_PATH $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download
wget http://downloads.sourceforge.net/project/picard/picard-tools/$VERSION/picard-tools-$VERSION.zip
unzip picard-tools-$VERSION.zip

# Move to install path
mv picard-tools-$VERSION $INSTALL_PATH
chmod -R g+w $INSTALL_PATH/picard-tools-$VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Java tool to manipulate BAMs \"
}
module-whatis \"MUGQIC - picard  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/picard/picard-tools-$VERSION
setenv          PICARD_HOME        \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/picard
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/picard/
rm -rf $INSTALL_DOWNLOAD
