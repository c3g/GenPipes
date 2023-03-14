#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ichorCNA
COMMIT=5bfc03e
VERSION=master-${COMMIT}
ARCHIVE=${SOFTWARE}-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/broadinstitute/${SOFTWARE}/archive/refs/heads/master.zip
# ARCHIVE_URL=https://github.com/broadinstitute/${SOFTWARE}/archive/refs/tags/v${VERSION}.tar.gz
SOFTWARE_DIR=${SOFTWARE}-${VERSION}
R_MODULE=mugqic/R_Bioconductor/4.2.2_3.16

build() {
  cd ${INSTALL_DOWNLOAD}

  rm -rf ${SOFTWARE}
  # git clone https://github.com/broadinstitute/${SOFTWARE}.git -b v${VERSION}
  git clone https://github.com/broadinstitute/${SOFTWARE}.git
  cd ${SOFTWARE}
  git checkout ${COMMIT}
  cd ..

  module load ${R_MODULE}
  mkdir -p ${INSTALL_DIR}/${SOFTWARE_DIR}

  cat << EOF > install_r_packages.R
install.packages("plyr", repos='http://cran.rstudio.org', lib="${INSTALL_DIR}/${SOFTWARE_DIR}")
BiocManager::install("HMMcopy", lib="${INSTALL_DIR}/${SOFTWARE_DIR}")
BiocManager::install("GenomeInfoDb", lib="${INSTALL_DIR}/${SOFTWARE_DIR}")
BiocManager::install("GenomicRanges", lib="${INSTALL_DIR}/${SOFTWARE_DIR}")
EOF

  Rscript --vanilla install_r_packages.R
  R CMD INSTALL -l ${INSTALL_DIR}/${SOFTWARE_DIR} ${SOFTWARE}
  cp -r ${SOFTWARE}/scripts ${INSTALL_DIR}/${SOFTWARE_DIR}/scripts
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

prereq mugqic/R_Bioconductor/4.0.3_3.12 mugqic/R_Bioconductor/4.1.0_3.13 mugqic/R_Bioconductor/4.2.1_3.15 mugqic/R_Bioconductor/4.2.2_3.16

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    R_LIBS              \$root
setenv          ICHORCNA_HOME       \$root/scripts
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@