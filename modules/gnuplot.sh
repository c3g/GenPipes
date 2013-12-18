MUGQIC_INSTALL_HOME="/sb/programs/analyste" # abacus

umask 0002
#cd $MUGQIC_INSTALL_HOME
##  Software dependencies installtion trace.


## Prepare
###################
################### gnuplot 4.6.4.
###################

## Install gnuplot itself
VERSION="4.6.4"
NAME=gnuplot

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$NAME/$NAME-$VERSION # where to install.
mkdir -p $MUGQIC_INSTALL_HOME/software/$NAME
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/bin

#wget http://sourceforge.net/projects/gnuplot/files/latest/download?source=dlp
#tar -xvf $NAME-$VERSION.tar.gz
#cd $NAME-$VERSION
#.configure --prefix=$MUGQIC_INSTALL_HOME/software/$NAME/$NAME-$VERSION
#make
#make install

# Module def file..
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Adds $NAME-$VERSION to your environment \"
}
module-whatis \"MUGQIC - Adds $NAME-$VERSION to your environment \"
                       
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/$NAME/$NAME-$VERSION
prepend-path    PATH               \$root/bin

" > $VERSION

## THEN--> Move module definition manually, and edit .version

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mkdir -p  $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$NAME
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$NAME

