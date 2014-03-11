#!/bin/sh

# rm -rf /mnt/parallel_scratch_mp2_wipe_on_august_2014/bourque/bourque_group/analyste/software/cd-hit
SOFTWARE="cd-hit"
VERSION="4.5.4-2011-03-07"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
wget http://www.bioinformatics.org/downloads/index.php/cd-hit/cd-hit-v$VERSION.tgz
tar xvf "cd-hit-v$VERSION.tgz"
cd "cd-hit-v$VERSION"
make -j8 openmp=yes


cd ..


# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .
mv -if cd-hit-v$VERSION.tgz $MUGQIC_INSTALL_HOME/archive/ 



# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ; 
}
module-whatis \"$SOFTWARE  \" ; 

set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/cd-hit-v$VERSION ;  
prepend-path    PATH                \$root ;
" > $VERSION


################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE


