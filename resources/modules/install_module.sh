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

download_archive() {
  INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
  mkdir -p $INSTALL_DOWNLOAD

  if [[ "$#" -eq 2 ]]
  then
      ARCHIVE_TMP=$2
      ARCHIVE_URL_PREFIX_TMP=$1
      ARCHIVE_URL_TMP=${ARCHIVE_URL_PREFIX_TMP}/${ARCHIVE_TMP}
  else
      ARCHIVE_TMP=$ARCHIVE
      ARCHIVE_URL_TMP=$ARCHIVE_URL
  fi

  # If archive was previously downloaded, use the local one, otherwise get it from remote site
  if [[ -f $ARCHIVE_DIR/$ARCHIVE_TMP ]]
  then
    echo "Archive $ARCHIVE_TMP already in $ARCHIVE_DIR/: using it..."
    cp -a $ARCHIVE_DIR/$ARCHIVE_TMP $INSTALL_DOWNLOAD/
  else
    echo "Archive $ARCHIVE_TMP not in $ARCHIVE_DIR/: downloading it..."
    wget --no-check-certificate $ARCHIVE_URL_TMP --output-document=$INSTALL_DOWNLOAD/$ARCHIVE_TMP
  fi
}

store_archive() {
  ARCHIVE_TMP=$1

  # Store archive if not already present
  if [[ ! -f $ARCHIVE_DIR/$ARCHIVE_TMP ]]
  then
    chmod -R ug+rwX,o+rX-w $INSTALL_DOWNLOAD/$ARCHIVE_TMP
    create_dir $ARCHIVE_DIR
    mv $INSTALL_DOWNLOAD/$ARCHIVE_TMP $ARCHIVE_DIR/
  fi
}

create_c3g_wrappers() {
  for i in `find $INSTALL_DIR/$SOFTWARE_DIR/ -type f -executable -exec file {} \; | grep ELF | grep -vP "\.so" | cut -d":" -f1`; do
    mv $i $i.raw;
    echo "$C3G_SYSTEM_LIBRARY/lib64/ld-linux-x86-64.so.2 --library-path $C3G_SYSTEM_LIBRARY/lib64:$C3G_SYSTEM_LIBRARY/lib64/mysql $i.raw \${@}" > $i;
  done
}

patch_c3g_binaries() {
  for i in `find $INSTALL_DIR/$SOFTWARE_DIR/ -type f -executable -exec file {} \; | grep ELF | grep -v "statically linked" | cut -d":" -f1`; do
    if readelf -l $i | grep go.build > /dev/null
    then
      echo "GO Done" > /dev/null
    elif [ ${i##*.} == "so" ] || [[ ${i##*/} =~ "so"*(\.[0-9])*$ ]]
    then
      $MUGQIC_INSTALL_HOME/software/PatchELF/PatchELF-0.9/bin/patchelf --set-rpath $C3G_SYSTEM_LIBRARY/usr/lib64/ $i
    else
      $MUGQIC_INSTALL_HOME/software/PatchELF/PatchELF-0.9/bin/patchelf --set-interpreter $C3G_SYSTEM_LIBRARY/lib64/ld-linux-x86-64.so.2 --set-rpath $C3G_SYSTEM_LIBRARY/usr/lib64/ $i
    fi
  done
}

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
elif [[ ${1:-} == MUGQIC_INSTALL_HOME_TMP ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME_TMP
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
fi

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE
ARCHIVE_DIR=${!INSTALL_HOME}/archive

# Set module directory path by lowercasing $INSTALL_HOME and removing '_install_home' in it
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME,,} | sed 's/_install_home//'`/$SOFTWARE

# Set path to C3G system libraries
C3G_SYSTEM_LIBRARY=/cvmfs/soft.mugqic/yum/centos7/1.0

echo "Installing $SOFTWARE version $VERSION in \$$INSTALL_HOME..."
echo

# Abort if software and/or module are already installed
if [[ -e $INSTALL_DIR/$SOFTWARE_DIR || -e $MODULE_DIR/$VERSION ]]
then
  echo "$INSTALL_DIR/$SOFTWARE_DIR and/or $MODULE_DIR/$VERSION already exist; please, delete them first!"
  exit 1
fi

create_dir $INSTALL_DIR
download_archive
build
patch_c3g_binaries

chmod -R ug+rwX,o+rX-w $INSTALL_DIR/$SOFTWARE_DIR

store_archive $ARCHIVE

# Deploy module
create_dir $MODULE_DIR
# Surround variable with "" since it contains a multiline text
module_file > $MODULE_DIR/$VERSION
# Default module version file
echo "create .version"
echo "\
#%Module1.0
set ModulesVersion \"$VERSION\"" > $MODULE_DIR/.version

echo "changing permission"
chmod ug+rwX,o+rX-w  $MODULE_DIR/$VERSION $MODULE_DIR/.version

# Clean up temporary installation files if any
cd
rm -rf $INSTALL_DOWNLOAD

echo
echo "$SOFTWARE version $VERSION has been successfully installed in \$$INSTALL_HOME"
if [[ $INSTALL_HOME == 'MUGQIC_INSTALL_HOME_DEV' ]]
then
  echo "To install module in production, type '$0 MUGQIC_INSTALL_HOME' (no '\$' before parameter)"
fi
