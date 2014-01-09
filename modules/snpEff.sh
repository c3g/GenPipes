
###################
################### SnpEff
###################
VERSION=3.4
# Replace "." in official version number by "_" in archive version number
ARCHIVE_VERSION=${VERSION//./_}
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/snpEff
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_PATH/archive $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download and extract
wget http://sourceforge.net/projects/snpeff/files/snpEff_v${ARCHIVE_VERSION}_core.zip
unzip snpEff_v${ARCHIVE_VERSION}_core.zip
mv snpEff $INSTALL_PATH/snpEff_$ARCHIVE_VERSION
cd $INSTALL_PATH
chmod -R g+w snpEff_$ARCHIVE_VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - SnpEff annotation tool \"
}
module-whatis \"MUGQIC - snpEff  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/snpEff/snpEff_$ARCHIVE_VERSION
setenv          SNPEFF_HOME        \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/snpEff
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/snpEff/
mv $INSTALL_DOWNLOAD/snpEff_v${ARCHIVE_VERSION}_core.zip $INSTALL_PATH/archive/
rm -rf $INSTALL_DOWNLOAD
