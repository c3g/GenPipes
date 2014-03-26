
###################
################### picard
###################
VERSION="1.1.4"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/mutect/
mkdir -p $INSTALL_PATH
# Download
cd $INSTALL_PATH
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/muTect-${VERSION}-bin.zip
unzip muTect-${VERSION}-bin.zip
# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Broads cancer snp caller \"
}
module-whatis \"MUGQIC - Broads cancer snp caller  \"
            
set             root              \$::env(MUGQIC_INSTALL_HOME)/software/mutect/muTect-${VERSION}
setenv          MUTECT_JAR        \$root/muTect-${VERSION}.jar
setenv          MUTECT_HOME       \$root/
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/mutect
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/mutect/




