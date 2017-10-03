#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=mirdeep2
VERSION=0.0.5
ARCHIVE=${SOFTWARE}-${VERSION}.zip
#ARCHIVE_URL=https://www.mdc-berlin.de/45995549/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/${SOFTWARE}_${VERSION//./_}.zip
ARCHIVE_URL=https://www.mdc-berlin.de/internet/38350089/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/older_miRDeep2_versions/mirdeep2_0_0_5.zip
SOFTWARE_DIR=$SOFTWARE-$VERSION

module_perl=mugqic/perl/5.22.1
module_bowtie=mugqic/bowtie/1.1.2

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
#  mv ${SOFTWARE} $SOFTWARE_DIR			# for VERSION=0.0.5
  mv ${SOFTWARE}_${VERSION//./_} $SOFTWARE_DIR	# for VERSION=0.0.8

  cd $SOFTWARE_DIR
  module load $module_perl $module_bowtie
  perl install.pl
#  for i in `ls *.pl`; do sed -ir 's/#\!\/usr\/bin\/perl/#\!\/usr\/bin\/env perl/' $i; done				# for VERSION=0.0.5
#  for i in `ls *.pl`; do sed -ir 's/#\!\/usr\/bin\/env perl -W/#\!\/usr\/bin\/env perl\nuse warnings;/' $i; done	# for VERSION=0.0.5
  for i in `ls bin/*.pl`; do sed -ir 's/#\!\/usr\/bin\/perl/#\!\/usr\/bin\/env perl/' $i; done				# for VERSION=0.0.8
  for i in `ls bin/*.pl`; do sed -ir 's/#\!\/usr\/bin\/env perl -W/#\!\/usr\/bin\/env perl\nuse warnings;/' $i; done	# for VERSION=0.0.8

  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
  rm -f $INSTALL_DIR/$SOFTWARE_DIR/bin/randfold
  rm -f $INSTALL_DIR/$SOFTWARE_DIR/bin/RNAfold
  mv $INSTALL_DIR/$SOFTWARE_DIR/essentials $INSTALL_DIR/    					# for VERSION=0.0.8
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
prepend-path    PATH                \$root/bin
prepend-path    PATH                \$root/../essentials/randfold-2.0
prepend-path    PATH                \$root/../essentials/ViennaRNA-1.8.4/install_dir/bin
module load $module_perl $module_bowtie
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
