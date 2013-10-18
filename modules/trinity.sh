
###################
################### Trinity
###################
VERSION=2013_08_14
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/trinity/
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_PATH/archive $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download and extract
wget http://sourceforge.net/projects/trinityrnaseq/files/trinityrnaseq_r$VERSION.tgz
tar xzvf trinityrnaseq_r$VERSION.tgz

# Compile
cd trinityrnaseq_r$VERSION
make -j8
cd ..

# Install
mv trinityrnaseq_r$VERSION $INSTALL_PATH
cd $INSTALL_PATH
chmod -R g+w trinityrnaseq_r$VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Trinity RNA assembler \"
}
module-whatis \"MUGQIC - Trinity  \"
            
set             root         \$::env(MUGQIC_INSTALL_HOME)/software/trinity/trinityrnaseq_r$VERSION
setenv          TRINITY_HOME \$root
prepend-path    PATH         \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trinity
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trinity/
mv $INSTALL_DOWNLOAD/trinityrnaseq_r$VERSION.tgz $INSTALL_PATH/archive/
rm -rf $INSTALL_DOWNLOAD
