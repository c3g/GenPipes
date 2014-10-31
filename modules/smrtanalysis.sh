#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# SMRT Analysis
#

SOFTWARE=smrtanalysis
VERSION_BASE=2.3.0.140936
VERSION_PATCH=p0
VERSION=$VERSION_BASE.$VERSION_PATCH

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV

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
ARCHIVE_BASE=${SOFTWARE}_$VERSION_BASE.run
ARCHIVE_PATCH=${SOFTWARE}-patch_$VERSION.run

for ARCHIVE in $ARCHIVE_BASE $ARCHIVE_PATCH
do
  # If archive was previously downloaded, use the local one, otherwise get it from remote site
  if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
  then
    echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
    cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
  else
    echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
    wget http://files.pacb.com/software/$SOFTWARE/${VERSION%\.*}/$ARCHIVE
  fi
done
# Extract and apply patch but NOT install SMRT Portal etc.
bash $ARCHIVE_BASE --extract-only --patchfile $ARCHIVE_PATCH

SOFTWARE_DIR=${SOFTWARE}_$VERSION

# Rename software default directory name with proper version number including patch
mv $SOFTWARE/install/${SOFTWARE}_$VERSION_BASE $SOFTWARE_DIR

# Default SGE cluster manager is not available; bash commands are run instead
sed -i 's/^CLUSTER_MANAGER = SGE/CLUSTER_MANAGER = BASH/' $SOFTWARE_DIR/analysis/etc/smrtpipe.rc

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX-w .
mv -i $SOFTWARE_DIR $INSTALL_DIR/
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
setenv          SEYMOUR_HOME        \$root
prepend-path    PATH                \$root/analysis/bin
puts            stderr              \"!!!===> Don't forget to source \\\${SEYMOUR_HOME}/etc/setup.sh  <===!!!\"
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
