MUGQIC_INSTALL_HOME="/sb/programs/analyste" # abacus

umask 0002
#cd $MUGQIC_INSTALL_HOME
##  Software dependencies installtion trace.
### Wed Dec 18 14:57:44 2013 


## Prepare
###################
################### gnuplot 4.6.4.
###################

## Install mem itself
VERSION="1.3"
NAME=memtime

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$NAME/$NAME-$VERSION # where to install.
mkdir -p $MUGQIC_INSTALL_HOME/software/$NAME
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/bin

wget http://www.update.uu.se/~johanb/memtime/memtime-1.3.tar.gz # Actually for our current version I pulled from the JGI repo from which I have access.
																# I believe there might have been some sligt modification to improve the code. Rob Egan is
																# maintaining the code on the JGI repo. rsegan@lbl.gov
tar -xvf $NAME-$VERSION.tar.gz
cd $NAME-$VERSION
.configure --prefix=$MUGQIC_INSTALL_HOME/software/$NAME/$NAME-$VERSION
make
make install

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
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$NAME/

