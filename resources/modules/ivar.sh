#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=ivar
VERSION=1.3
ARCHIVE=$SOFTWARE-${VERSION}.tar.gz
ARCHIVE_URL=https://github.com/andersen-lab/${SOFTWARE}/archive/v${VERSION}.tar.gz
SOFTWARE_DIR=$SOFTWARE-$VERSION

#module () {   
#    eval `/usr/bin/modulecmd sh $*`
#}

build() {
  cd $INSTALL_DOWNLOAD
  tar -xzvf $ARCHIVE

  cd $SOFTWARE_DIR
  module load mugqic/htslib/1.11
  ./autogen.sh
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR --with-hts=$HTSLIB_HOME
  make
  make install
  if $MUGQIC_INSTALL_HOME/software/patchelf/patchelf-0.9/bin/patchelf --print-rpath $INSTALL_DIR/$SOFTWARE_DIR/bin/ivar
  then
    $MUGQIC_INSTALL_HOME/software/patchelf/patchelf-0.9/bin/patchelf --set-rpath $($MUGQIC_INSTALL_HOME/software/patchelf/patchelf-0.9/bin/patchelf --print-rpath $INSTALL_DIR/$SOFTWARE_DIR/bin/ivar):$HTSLIB_LIBRARY_DIR $INSTALL_DIR/$SOFTWARE_DIR/bin/ivar
  else
    $MUGQIC_INSTALL_HOME/software/patchelf/patchelf-0.9/bin/patchelf --set-rpath $HTSLIB_LIBRARY_DIR $INSTALL_DIR/$SOFTWARE_DIR/bin/ivar
  fi

  cd bin
  wget https://raw.githubusercontent.com/c3g/CoVSeQ/dev/ivar_variants_to_vcf.py
  chmod a+x ivar_variants_to_vcf.py
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
if { [ module-info mode load ] } {
    puts stderr \"WARNING : iVar needs both samtools and htslib to be accessible in the envoronement.\"
    puts stderr \"We recommend loading mugqic/samtools/1.11 and mugqic/htslib/1.11 \"
}
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
<<<<<<< HEAD

=======
>>>>>>> 0dbd8958adc117d03f7623c692579e981b4e0029
