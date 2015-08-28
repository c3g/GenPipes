#!/bin/sh

SOFTWARE="popoolation2"
VERSION="1201"  ## Checked out svn revision 189
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:+
wget  http://downloads.sourceforge.net/project/popoolation2/popoolation2_$VERSION.zip -O popoolation2_$VERSION.zip
unzip popoolation2_$VERSION.zip

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
chmod -R o+rX .
mv -i $SOFTWARE"_"$VERSION $INSTALL_PATH/$SOFTWARE-$VERSION


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  
}
module-whatis \"$SOFTWARE  is designed to analyze next generation sequencing data of pooled genomic DNA.\" ;  

prereq  mugqic/perl/5.18.2                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$SOFTWARE-$VERSION ;  
prepend-path    PATH                \$root;  
prepend-path    PERL5LIB            \$root/Modules; 
setenv          POPOOLATION_PATH    \$root;
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD


