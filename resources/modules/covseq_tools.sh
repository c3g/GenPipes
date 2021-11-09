#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=CoVSeQ_tools
VERSION=1.2.0
ARCHIVE=${SOFTWARE,,}-$VERSION.tar.gz
ARCHIVE_URL=https://bitbucket.org/mugqic/${SOFTWARE}/downloads/$ARCHIVE
SOFTWARE_DIR=${SOFTWARE,,}-$VERSION

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
prepend-path    PATH                   \$root/run_reporting
prepend-path    PATH                   \$root/illumina_metrics
prepend-path    PATH                   \$root/full_reporting
setenv          RUN_REPORT             \$root/run_reporting/run_report.Rmd
setenv          RUN_REPORT_FREEBAYES   \$root/run_reporting/run_report_freebayes.Rmd
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
