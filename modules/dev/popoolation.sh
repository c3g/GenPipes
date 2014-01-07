#!/bin/sh

SOFTWARE="popoolation"
VERSION="0.218"  ## Checked out svn revision 218
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
svn checkout http://popoolation.googlecode.com/svn $SOFTWARE

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
chmod -R o+rX .
mv -i $SOFTWARE-$VERSION $INSTALL_PATH  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION IF DIFFERENT
zip -r $SOFTWARE-$VERSION.zip popoolation/*
mv -i $SOFTWARE-$VERSION.zip $MUGQIC_INSTALL_HOME_DEV/archive  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2) IF DIFFERENT

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  
}
module-whatis \"$SOFTWARE  is designed to analyze next generation sequencing data of pooled genomic DNA.\" ;  
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;  
prepend-path    PATH                \$root/popoolation/ ;  
prepend-path    PATH                \$root/other_tools/bin ; 
setenv          PERL5LIB            \$root/popoolation/Modules ; 
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
