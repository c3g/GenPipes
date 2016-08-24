#!/bin/bash

#
# ea-utils
#

SOFTWARE=ea-utils
VERSION=1.1.2-537
INSTALL_PATH=$MUGQIC_INSTALL_HOME_TMP/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget http://ea-utils.googlecode.com/files/$SOFTWARE.$VERSION.tar.gz
tar zxvf $SOFTWARE.$VERSION.tar.gz
cd $SOFTWARE.$VERSION
# Guillimin requires GSL module and g++ version 4.7 to compile
if [[ `hostname` == lg-* || `dnsdomainname` == guillimin.clumeq.ca ]]
then
  source /etc/profile.d/modules.sh
  module load gcc/4.7.2 GSL/1.15
fi
PREFIX=$INSTALL_PATH/$SOFTWARE.$VERSION make install

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX $INSTALL_PATH/$SOFTWARE.$VERSION .
mv -i $SOFTWARE.$VERSION.tar.gz $MUGQIC_INSTALL_HOME_TMP/archive

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE  \"

set             root                \$::env(MUGQIC_INSTALL_HOME_TMP)/software/$SOFTWARE/$SOFTWARE.$VERSION
prepend-path    PATH                \$root/bin
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_TMP/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_TMP/modulefiles/mugqic/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
