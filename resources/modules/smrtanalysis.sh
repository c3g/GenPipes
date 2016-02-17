#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=smrtanalysis
VERSION_BASE=2.3.0.140936
VERSION_PATCH=p5
VERSION=$VERSION_BASE.$VERSION_PATCH
ARCHIVE=${SOFTWARE}_$VERSION_BASE.run
ARCHIVE_PATCH=${SOFTWARE}-patch_$VERSION.run
ARCHIVE_URL_PREFIX=http://files.pacb.com/software/$SOFTWARE/${VERSION_BASE%\.*}
ARCHIVE_URL=$ARCHIVE_URL_PREFIX/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE}_$VERSION
#https://files.pacb.com/software/smrtanalysis/2.3.0/smrtanalysis-patch_2.3.0.140936.p4.run

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {

  cd $INSTALL_DOWNLOAD

  # Download main source
  download_archive
  # Download patch too
  download_archive $ARCHIVE_URL_PREFIX $ARCHIVE_PATCH 

  # Bash cannot run patch if not executable
  chmod +x $ARCHIVE $ARCHIVE_PATCH
  # Extract and apply patch but NOT install SMRT Portal etc.
  bash $ARCHIVE --extract-only --patchfile $ARCHIVE_PATCH

  # Rename software default directory name with proper version number including patch
  mv $SOFTWARE/install/${SOFTWARE}_$VERSION_BASE $SOFTWARE_DIR

  # Default SGE cluster manager is not available; bash commands are run instead
  sed -i 's/^CLUSTER_MANAGER = SGE/CLUSTER_MANAGER = BASH/' $SOFTWARE_DIR/analysis/etc/smrtpipe.rc

  # Install software
  mv -i $SOFTWARE_DIR $INSTALL_DIR/

  # Archive patch file as well
  store_archive $ARCHIVE_PATCH
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
setenv          SEYMOUR_HOME        \$root
prepend-path    PATH                \$root/analysis/bin
puts            stderr              \"!!!===> Don't forget to source \\\${SEYMOUR_HOME}/etc/setup.sh  <===!!!\"
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
