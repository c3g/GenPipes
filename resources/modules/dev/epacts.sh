#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# NOTES:
# - The script 	$INSTALL_DIR/$SOFTWARE_DIR/bin/epacts download downloads ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz with perl NET::FTP, this fails on abacus for some reason
# - http://csg.sph.umich.edu/kang/epacts/download/$ARCHIVE was originally blocketd on abacus, now fixed.


# (See http://genome.sph.umich.edu/wiki/EPACTS for comprehensive documentation)
#
# == EPACTS Installation Details  ==
#
# If you want to use EPACTS in an Ubuntu platform, following the step below
#
# * Download EPACTS source distribution at http://www.sph.umich.edu/csg/kang/epacts/download/EPACTS-3.0.0.tar.gz (99MB)
# * Uncompress EPACTS package, and install the package using the following set of commands
#
#   tar xzvf EPACTS-$(VERSION).tar.gz
#   cd EPACTS-$(VERSION)
#   ./configure --prefix=/path/to/install
#   (If you have libraries in non-standard directory, please try
#    ./configure --prefix=/path/to/install LDFLAGS=-L/path/to/library to include the directory containing libR.so)
#   make
#   make install
#
# (Important Note: '''make sure to specify --prefix=/path/to/install''' to avoid installing to the default path /usr/local/, which you may not have the permission. /home/your_userid/epacts might be a good one, if you are not sure where to install)
#
# * Now ${EPACTS_DIR} represents the '/path/to/install' directory
#
# * Download the reference FASTA files from 1000 Genomes FTP automatically by running the following commands
#
#   ${EPACTS_DIR}/bin/epacts download
#
#  (For advanced users, to save time for downloading the FASTA files (~900MB), you may copy a local copy of GRCh37 FASTA file and the index file to ${EPACTS_DIR}/share/EPACTS/)
#
# *Perform a test run by running the following command
#
#   ${EPACTS_DIR}/bin/test_run_epacts.sh


PREREQS="mugqic/gnuplot/4.6.6 mugqic/R_Bioconductor/3.1.2_3.0"
SOFTWARE="EPACTS"
VERSION="3.2.6"  
ARCHIVE="$SOFTWARE-$VERSION.tar.gz"  
ARCHIVE_URL="http://csg.sph.umich.edu/kang/epacts/download/$ARCHIVE"
SOFTWARE_DIR=$SOFTWARE-$VERSION 

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE 
  cd $SOFTWARE_DIR
	module load $PREREQS
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR  
  make 
  make install  
	
	$INSTALL_DIR/$SOFTWARE_DIR/bin/epacts download
	$INSTALL_DIR/$SOFTWARE_DIR/bin/test_run_epacts.sh &> $INSTALL_DIR/$SOFTWARE_DIR/test_run_epacts.sh.out   # $HOME/test/bin/test_run_epacts.sh
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

prereq mugqic/gnuplot/4.6.6 mugqic/R_Bioconductor/3.1.2_3.0

set             root                $INSTALL_DIR/$SOFTWARE_DIR ;
prepend-path    PATH                \$root/bin ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
