#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# MACS
#

SOFTWARE=MACS
VERSION=2.0.10.09132012
#VERSION=1.4.2
PYTHON_VERSION=2.7.8

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX-w $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
if [[ $VERSION == "1.4.2" ]]
then
  ARCHIVE=$SOFTWARE-$VERSION-1.tar.gz
  URL=https://github.com/downloads/taoliu/MACS/$ARCHIVE
  SOFTWARE_DIR=$SOFTWARE-$VERSION
elif [[ $VERSION == "2.0.10.09132012" ]]
then
  ARCHIVE=${SOFTWARE}2-$VERSION.tar.gz
  URL=https://pypi.python.org/packages/source/M/MACS2/$ARCHIVE
  SOFTWARE_DIR=${SOFTWARE}2-$VERSION
fi

# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget --no-check-certificate $URL
fi
tar zxvf $ARCHIVE

cd $SOFTWARE_DIR
module load mugqic/python/$PYTHON_VERSION
python setup.py install --prefix $INSTALL_DIR/$SOFTWARE_DIR

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX-w . $INSTALL_DIR/$SOFTWARE_DIR
# Store archive if not already present or if different from the previous one
if [[ ! -f ${!INSTALL_HOME}/archive/$ARCHIVE || `diff ${!INSTALL_HOME}/archive/$ARCHIVE $ARCHIVE` ]]
then
  mv -i $ARCHIVE ${!INSTALL_HOME}/archive/
fi

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    PYTHONPATH         \$root/lib/python2.7/site-packages
setenv          MACS_BIN            \$root/bin
setenv          MACS_LIB            \$root/lib/python2.7/site-packages
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by lowercasing $INSTALL_HOME and removing '_install_home' in it
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME,,} | sed 's/_install_home//'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX-w $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX-w $VERSION .version
mv $VERSION .version $MODULE_DIR/

# Clean up temporary installation files if any
cd
rm -rf $INSTALL_DOWNLOAD
