#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=VarDictJava 
VERSION=1.4.9
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/AstraZeneca-NGS/$SOFTWARE/archive/v$VERSION.tar.gz  
SOFTWARE_DIR=$SOFTWARE-$VERSION 

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE  

  #Unzip compile version
  cd $SOFTWARE_DIR/dist
  if [[ $VERSION == "1.4.10" ]]
  then
    wget https://github.com/AstraZeneca-NGS/$SOFTWARE/files/758679/VarDict-1.4.10.zip
  elif [[ $VERSION == 1.4.9 ]]
  then
    wget https://github.com/AstraZeneca-NGS/$SOFTWARE/files/748361/VarDict-1.4.9.zip
  fi

  unzip VarDict-${VERSION}.zip
  mv VarDict-${VERSION}/* ../
  
  cd ../
  
  #Install Vardict perl support files  
  rm -R VarDict
  git clone https://github.com/AstraZeneca-NGS/VarDict.git
 
  # Install software
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
prepend-path    PATH                \$root/bin ; 
setenv          VARDICT_HOME        \$root ;
setenv          VARDICT_BIN         \$root/VarDict ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
