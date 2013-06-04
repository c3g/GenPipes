
###################
################### picard
###################
VERSION="2013-02-25";
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/trinity/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download
wget http://sourceforge.net/projects/trinityrnaseq/files/trinityrnaseq_r${VERSION}.tgz
tar xzvf trinityrnaseq_r${VERSION}.tgz
cd trinityrnaseq_r${VERSION}
make
# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Trinity RNA assembler \"
}
module-whatis \"MUGQIC - Trinity  \"
            
set             root      \$::env(MUGQIC_INSTALL_HOME)/software/trinity/trinityrnaseq_r${VERSION}
prepend-path    PATH      \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trinity
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trinity/




