
###################
################### picard
###################
VERSION="3.2"
VERSION_DL="3_2"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/snpEff/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download
wget http://sourceforge.net/projects/snpeff/files/snpEff_v${VERSION_DL}_core.zip
unzip snpEff_v${VERSION_DL}_core.zip
mv snpEff snpEff_${VERSION_DL}
# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - SnpEff annotation tool \"
}
module-whatis \"MUGQIC - snpEff  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/snpEff/snpEff_${VERSION_DL}
setenv          SNPEFF_HOME        \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/snpEff
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/snpEff/




