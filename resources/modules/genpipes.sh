#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=genpipes
VERSION=4.1.2
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://bitbucket.org/mugqic/genpipes/downloads/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

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
prepend-path    PATH                  \$root/pipelines/ampliconseq
prepend-path    PATH                  \$root/pipelines/covseq
prepend-path    PATH                  \$root/pipelines/chipseq
prepend-path    PATH                  \$root/pipelines/dnaseq
prepend-path    PATH                  \$root/pipelines/dnaseq_high_coverage
prepend-path    PATH                  \$root/pipelines/epiqc
prepend-path    PATH                  \$root/pipelines/hicseq
prepend-path    PATH                  \$root/pipelines/methylseq
prepend-path    PATH                  \$root/pipelines/nanopore
prepend-path    PATH                  \$root/pipelines/nanopore_covseq
prepend-path    PATH                  \$root/pipelines/rnaseq
prepend-path    PATH                  \$root/pipelines/rnaseq_denovo_assembly
prepend-path    PATH                  \$root/pipelines/rnaseq_light
prepend-path    PATH                  \$root/pipelines/hicseq
prepend-path    PATH                  \$root/pipelines/tumor_pair
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
