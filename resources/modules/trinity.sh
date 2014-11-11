#!/bin/bash

#
# Trinity
#

SOFTWARE=trinity
VERSION=20140413p1

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
ARCHIVE=${SOFTWARE}rnaseq_r$VERSION.tar.gz
wget http://sourceforge.net/projects/trinityrnaseq/files/$ARCHIVE
tar zxvf $ARCHIVE

SOFTWARE_DIR=${SOFTWARE}rnaseq_r$VERSION
cd $SOFTWARE_DIR
# Set Makefile flags to compile jellyfish statically
sed -i 's/cd trinity-plugins\/jellyfish \&\& \.\/configure.*/cd trinity-plugins\/jellyfish \&\& .\/configure CC=gcc CXX=g++ --enable-static --disable-shared --prefix=`pwd` \&\& $(MAKE) LDFLAGS=-all-static AM_CPPFLAGS="-Wall -Wnon-virtual-dtor -I"`pwd`/' Makefile
# Disable "-" to "_" substitution in abundance matrix since we don't want sample names to be modified
sed -ri 's/(\$file =\~ s\/\\-\/_\/g;)/#\1/' util/abundance_estimates_to_matrix.pl
make

# Add permissions and install archive
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX .
mv -i $SOFTWARE_DIR $INSTALL_DIR
mv -i $ARCHIVE ${!INSTALL_HOME}/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE\" ;

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
setenv          TRINITY_HOME        \$root ;
prepend-path    PATH                \$root ;
prepend-path    PATH                \$root/util ;
prepend-path    PATH                \$root/util/RSEM_util ;
prepend-path    PATH                \$root/Analysis/DifferentialExpression ;
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by removing '_INSTALL_HOME' in $INSTALL_HOME and lowercasing the result
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME/_INSTALL_HOME/} | tr '[:upper:]' '[:lower:]'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX $VERSION .version
mv $VERSION .version $MODULE_DIR

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
