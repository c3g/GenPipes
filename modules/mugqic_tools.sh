#!/bin/bash

#
# MUGQIC Tools
#

SOFTWARE=mugqic_tools
VERSION=1.9
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget https://bitbucket.org/mugqic/$SOFTWARE/get/$VERSION.tar.gz
# Rename archive
mv $VERSION.tar.gz $SOFTWARE-$VERSION.tar.gz
tar zxvf $SOFTWARE-$VERSION.tar.gz
# Rename mugqic_tools directory since original bitbucket name contains the commit number instead of version
mv mugqic-mugqic_tools* $SOFTWARE-$VERSION

# Add permissions and install software
cd $INSTALL_DOWNLOAD 
chmod -R ug+rwX .
chmod -R o+rX .
mv -i $SOFTWARE-$VERSION $INSTALL_PATH
mv -i $SOFTWARE-$VERSION.tar.gz $MUGQIC_INSTALL_HOME/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;
prepend-path    PATH                \$root/tools ;
prepend-path    PATH                \$root/perl-tools ;
prepend-path    PATH                \$root/R-tools ;
prepend-path    PATH                \$root/python-tools ;
prepend-path    PATH                \$root/RRNATagger-tools ;
prepend-path    PERL5LIB            \$root/perl-tools ;
setenv          R_TOOLS             \$root/R-tools ;
setenv          PERL_TOOLS          \$root/perl-tools ;
setenv          PYTHON_TOOLS        \$root/python-tools ;
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Well... here, module directory is named "tools" instead of "mugqic_tools" for aesthetical reasons

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
MODULE_DIR=$MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools
mkdir -p $MODULE_DIR
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MODULE_DIR

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
