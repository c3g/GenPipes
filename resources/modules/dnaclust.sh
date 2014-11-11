#!/bin/bash

MUGQIC_INSTALL_HOME="/sb/programs/analyste" # abacus
# MUGQIC_OPT="/mnt/scratch_mp2/bourque/bourque_group/software/areas/genomics/opt/" # Mammouth (NO GOOD for now)
#MUGQIC_INSTALL_HOME="/software/areas/genomics" # Guillimin


umask 0002
#cd $MUGQIC_INSTALL_HOME
##  Software dependencies installtion trace.
### Mon 14 Jan 20:21:44 2013 


## Prepare
###################
################### dnaclust parallel release 3.
###################

## Install dnaclust itself
VERSION="3"
NAME=dnaclust
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$NAME/$NAME-$VERSION # where to install.
mkdir -p $MUGQIC_INSTALL_HOME/software/$NAME
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/bin
wget http://sourceforge.net/projects/dnaclust/files/parallel_release_$VERSION/dnaclust_linux_release$VERSION.zip/download
unzip dnaclust_linux_release$VERSION.zip
cd ${NAME}_linux_release$VERSION
cp * $INSTALL_PATH/bin/

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

