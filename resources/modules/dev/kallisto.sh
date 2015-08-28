#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE="kallisto" 
VERSION="0.42.2.1" 
ARCHIVE="v$VERSION.tar.gz"
ARCHIVE_URL="https://github.com/pachterlab/$SOFTWARE/archive/v$VERSION.tar.gz"
SOFTWARE_DIR=$SOFTWARE-$VERSION  


# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE  ## TO BE MODIFIED WITH SPECIFIC COMMAND

	# cd ..; rm -rf build; mkdir build; cd build

  cd $SOFTWARE_DIR
  mkdir build
 	cd build
	
	module load  itgenome/gcc/4.8.5  # abacus
	
	# module load cmake/2.8.8 hdf5/1.8.9  # mammouth
	# - CMake version >= 2.8.12
	#     - Can be installed via homebrew: `brew install cmake`
	# - zlib (should be installed on OSX >= 10.9)
	# - HDF5 C libraries
	#
	
  cmake -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_DIR/$SOFTWARE_DIR ..
	make
	make install
	
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

module load itgenome/gcc/4.8.5
set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin ; 
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
