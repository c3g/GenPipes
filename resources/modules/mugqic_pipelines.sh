#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=mugqic_pipelines
VERSION=2.2.0
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://bitbucket.org/mugqic/$SOFTWARE/downloads/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE

  # Install software
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                   $INSTALL_DIR/$SOFTWARE_DIR
setenv          MUGQIC_PIPELINES_HOME \$root
prepend-path    PATH                  \$root/utils
prepend-path    PATH                  \$root/pipelines/chipseq
prepend-path    PATH                  \$root/pipelines/dnaseq
prepend-path    PATH                  \$root/pipelines/dnaseq_high_coverage
prepend-path    PATH                  \$root/pipelines/illumina_run_processing
prepend-path    PATH                  \$root/pipelines/pacbio_assembly
prepend-path    PATH                  \$root/pipelines/rnaseq
prepend-path    PATH                  \$root/pipelines/rnaseq_denovo_assembly
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
