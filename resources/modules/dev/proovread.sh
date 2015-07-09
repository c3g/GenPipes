#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

SOFTWARE="proovread" 
VERSION="2.12"
ARCHIVE="$SOFTWARE-$VERSION.tar.gz"
ARCHIVE_URL="https://github.com/BioInf-Wuerzburg/proovread/archive/$SOFTWARE-$VERSION.tar.gz"
SOFTWARE_DIR=$SOFTWARE-$VERSION 

module load mugqic/perl/5.18.2
cpan -i Log::Log4perl


build() {
  cd $INSTALL_DOWNLOAD

	git clone --recursive https://github.com/BioInf-Wuerzburg/proovread
	mv $SOFTWARE $SOFTWARE_DIR
	cd $SOFTWARE_DIR/util/bwa
	make -j8

  # Install software
  cd $INSTALL_DOWNLOAD  ## TO BE ADDED AND MODIFIED IF NECESSARY
  mv -i $SOFTWARE_DIR $INSTALL_DIR/  ## TO BE ADDED AND MODIFIED IF NECESSARY
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"
module load mugqic/perl/5.18.2 mugqic_dev/samtools/1.2 mugqic/blast/2.2.29+
prepend-path PERL5LIB /software/areas/genomics/perl5libs/lib/site_perl/5.18.2/
set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
