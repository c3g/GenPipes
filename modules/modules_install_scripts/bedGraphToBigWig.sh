

###################
################### bedGraphToBigWig
###################
VERSION="v4"
NAME=bedGraphToBigWig # same could apply to all ucsc tools
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$NAME/$NAME-$VERSION
mkdir -p $INSTALL_PATH
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/$NAME
chmod +x $NAME
mv $NAME $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $NAME \"
}
module-whatis \"MUGQIC - $NAME \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/$NAME/$NAME-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$NAME
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$NAME/




