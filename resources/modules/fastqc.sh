#!/bin/bash

###################
################### FASTQC
###################
VERSION="0.11.2"
if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
  MODULE=mugqic
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
  MODULE=mugqic_dev
fi

INSTALL_PATH=${!INSTALL_HOME}/software/fastqc/fastqc_v"$VERSION" # where to install...

mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v"$VERSION".zip
unzip fastqc_v"$VERSION".zip
rm fastqc_v"$VERSION".zip
chmod +x FastQC/fastqc
chmod -R g+w $INSTALL_PATH
sed -i s,"#\!/usr/bin/perl.*,#\!/usr/bin/env perl,g" FastQC/fastqc

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - FASTQC \"
}
module-whatis \"MUGQIC -FASTQC \"
                      
set             root                $INSTALL_PATH/FastQC
prepend-path    PATH                \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p ${!INSTALL_HOME}/modulefiles/$MODULE/fastqc
mv .version $VERSION ${!INSTALL_HOME}/modulefiles/$MODULE/fastqc/


