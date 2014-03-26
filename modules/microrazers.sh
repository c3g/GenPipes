#!/bin/bash
SOFTWARE=seqan  
VERSION=1.4.1
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget http://packages.seqan.de/seqan-apps/seqan-apps-$VERSION-Linux-x86_64.tar.bz2
tar -xvf seqan-apps-$VERSION-Linux-x86_64.tar.bz2

# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .
mv -i seqan-apps-$VERSION-Linux-x86_64.tar.bz2 $MUGQIC_INSTALL_HOME/archive 



# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  ## TO BE MODIFIED WITH DETAILED DESCRIPTION IF ANY
}
module-whatis \"$SOFTWARE  \" ;  ## TO BE MODIFIED WITH DETAILED DESCRIPTION IF ANY
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/seqan-apps-$VERSION-Linux-x86_64 ;  
prepend-path    PATH                \$root/bin ;  
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


