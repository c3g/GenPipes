#!/bin/sh

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# You should probably also delete this commented-out header and the ## comments
################################################################################


#
# Software_name  ## TO BE MODIFIED WITH e.g. BLAST, HMMER, SAMtools, etc.
#

SOFTWARE=tigra  ## TO BE MODIFIED WITH e.g. blast, hmmer, samtools, etc.
VERSION=0.1  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE/$VERSION
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
git clone --recursive https://github.com/genome/tigra-sv.git

cd tigra-sv
mkdir build
cd build

wget "http://downloads.sourceforge.net/project/samtools/samtools/0.1.6/samtools-0.1.6.tar.bz2"
tar -jxvf samtools-0.1.6.tar.bz2
cd samtools-0.1.6
make
export SAMTOOLS_ROOT=$(pwd)
cd ..
cp samtools-0.1.6/* ../src/exe/tigra-sv/

module add mugqic_dev/cmake/2.8.12.1-Linux-i386

#TODO edit ../CMakeLists.txt
#set(SAMTOOLS_INC /lb/project/mugqic/analyste_dev/software/tigra/0.1/tigra-sv/build/samtools-0.1.6)
#set(SAMTOOLS_LIB /lb/project/mugqic/analyste_dev/software/tigra/0.1/tigra-sv/build/samtools-0.1.6)
cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH

make

make install

cd ..

# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  
}
module-whatis \"$SOFTWARE  \" ;  
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$VERSION/ ;  
prepend-path    PATH                \$root/tigra-sv/build/bin ;  
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE

# Clean up temporary installation files if any
#rm -rf $INSTALL_DOWNLOAD

