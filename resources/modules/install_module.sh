#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

create_dir() {
  DIR=$1

  # Create directory with permissions if necessary
  if [[ ! -d $DIR ]]
  then
    mkdir -p $DIR
    chmod ug+rwX,o+rX-w $DIR
  fi
}

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
fi

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE
ARCHIVE_DIR=${!INSTALL_HOME}/archive

# Set module directory path by lowercasing $INSTALL_HOME and removing '_install_home' in it
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME,,} | sed 's/_install_home//'`/$SOFTWARE

echo "Installing $SOFTWARE version $VERSION in \$$INSTALL_HOME..."
echo

# Abort if software and/or module are already installed
if [[ -e $INSTALL_DIR/$SOFTWARE_DIR || -e $MODULE_DIR/$VERSION ]]
then
  echo "$INSTALL_DIR/$SOFTWARE_DIR and/or $MODULE_DIR/$VERSION already exist; please, delete them first!"
  exit 1
fi

create_dir $INSTALL_DIR

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir -p $INSTALL_DOWNLOAD

# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f $ARCHIVE_DIR/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in $ARCHIVE_DIR/: using it..."
  cp -a $ARCHIVE_DIR/$ARCHIVE $INSTALL_DOWNLOAD/
else
  echo "Archive $ARCHIVE not in $ARCHIVE_DIR/: downloading it..."
  wget --no-check-certificate $ARCHIVE_URL --output-document=$INSTALL_DOWNLOAD/$ARCHIVE
fi

build $ARCHIVE

chmod -R ug+rwX,o+rX-w $INSTALL_DIR/$SOFTWARE_DIR $INSTALL_DOWNLOAD/$ARCHIVE

# Store archive if not already present
if [[ ! -f $ARCHIVE_DIR/$ARCHIVE ]]
then
  create_dir $ARCHIVE_DIR
  mv $INSTALL_DOWNLOAD/$ARCHIVE $ARCHIVE_DIR/
fi

# Deploy module
create_dir $MODULE_DIR
# Surround variable with "" since it contains a multiline text
module_file > $MODULE_DIR/$VERSION
# Default module version file
echo "\
#%Module1.0
set ModulesVersion \"$VERSION\"" > $MODULE_DIR/.version

chmod ug+rwX,o+rX-w $MODULE_DIR/$VERSION $MODULE_DIR/.version

# Clean up temporary installation files if any
cd
rm -rf $INSTALL_DOWNLOAD

echo
echo "$SOFTWARE version $VERSION has been successfully installed in \$$INSTALL_HOME"
if [[ $INSTALL_HOME == 'MUGQIC_INSTALL_HOME_DEV' ]]
then
  echo "To install module in production, type '$0 MUGQIC_INSTALL_HOME' (no '\$' before parameter)"
fi
