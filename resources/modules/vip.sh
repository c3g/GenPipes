#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=VIP
VERSION=1.0.1
ARCHIVE=${SOFTWARE}-$VERSION.zip
ARCHIVE_URL=https://github.com/keylabivdc/${SOFTWARE}/archive/master.zip
SOFTWARE_DIR=${SOFTWARE}-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  mkdir -p $SOFTWARE_DIR/bin
  mkdir -p $SOFTWARE_DIR/data

  cd $SOFTWARE_DIR/bin
  git clone --recursive git://github.com/keylabivdc/${SOFTWARE}.git
  cd VIP
  chmod 755 *
  chmod 755 *
  find -maxdepth 1 -type f -exec cp {} $INSTALL_DOWNLOAD/$SOFTWARE_DIR/bin/ \;
  mv edirect $INSTALL_DOWNLOAD/$SOFTWARE_DIR/
  mv installer $INSTALL_DOWNLOAD/$SOFTWARE_DIR/
  cd $INSTALL_DOWNLOAD/$SOFTWARE_DIR/bin
  rm -rf VIP
  PATH=$PATH:$INSTALL_DOWNLOAD/$SOFTWARE_DIR/edirect
 
  cd $INSTALL_DOWNLOAD/$SOFTWARE_DIR/ 
  sed -i "s|sudo||g" installer/db_installer.sh
  sed -i "s|/usr/VIP|$INSTALL_DIR/$SOFTWARE_DIR|g" installer/db_installer.sh
  sed -i "s|gzip -d -k|gzip -d|g" installer/db_installer.sh
  sed -i "s|gzip -dc -k|gzip -dc|g" installer/db_installer.sh
  bash installer/db_installer.sh -r $INSTALL_DOWNLOAD/$SOFTWARE_DIR/data

  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
prepend-path    PATH                \$root/bin
prepend-path    PATH                \$root/edirect
prepend-path    PATH                \$root/data
module load mugqic/bowtie2/2.2.9 mugqic/velvet/1.2.10 mugqic/mafft/7.310 mugqic/picard/2.0.1 mugqic/python/2.7.13 mugqic/perl/5.22.1
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@

