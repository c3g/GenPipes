#!/bin/sh

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


#
# Perl
#

SOFTWARE=perl   
VERSION=5.18.2  
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget http://www.cpan.org/src/5.0/$SOFTWARE-$VERSION.tar.gz     
tar -xvf $SOFTWARE-$VERSION.tar.gz                             
cd $SOFTWARE-$VERSION                                          
./Configure -des -Dusethreads -Dprefix=$INSTALL_PATH/$SOFTWARE-$VERSION     
make  
make install

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
chmod -R o+rX .
mv -i $SOFTWARE-$VERSION.tar.gz $MUGQIC_INSTALL_HOME/archive 

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE-$VERSION \" ;  
}
module-whatis \"$SOFTWARE-$VERSION  \" ;  
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION
prepend-path    PATH                \$root/bin                                                           
prepend-path    PERL5LIB            \$root/lib
prepend-path    PERL5LIB            \$root/lib/site_perl/5.18.2
prepend-path    PERL5LIB            \$root/lib/5.18.2
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

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
