#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=Basecall_T7_GPU
VERSION=1.4.20.71
ARCHIVE=${SOFTWARE}_${VERSION}_CentOS.run
# Empty URL triggers the usage of the installer in archive since it's not
# publicly available
ARCHIVE_URL=
# Original installer must be copied to the archive directory:
# $MUGQIC_INSTALL_HOME/archive/
T7_BASECALL_ORIGINAL=/nb/Research/MGISeq/T7/R1100600200054/Basecall_T7_GPU_1.4.20.71_CentOS.run
SOFTWARE_DIR=$VERSION

# TODO make this script copy the installer in the archive directory, before
# variable statements

build() {
    cd $INSTALL_DIR
    # Adapt installer, to sudoless run and replace target installPath
    cp $ARCHIVE_DIR/$ARCHIVE_TMP $ARCHIVE.tmp
    sed -i "{s|sudo ||g; s|installPath=/opt/basecall|installPath=$INSTALL_DIR|; s|runningPath=/home/basecall/proc|runningPath=$INSTALL_DIR|}" $ARCHIVE.tmp
    chmod +x $ARCHIVE.tmp
    # Install software
    ./$ARCHIVE.tmp
    mv ${ARCHIVE%.*} $VERSION
    rm $ARCHIVE.tmp
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
    puts stderr \"\tMUGQIC - $SOFTWARE\n\tThe basecaller version in use must\
    match the version specified in the flagfile of the flowcell or run.\"
}
module-whatis \"$SOFTWARE\"
if { [ module-info mode load ] } {
    puts stderr \"Info: The T7 basecaller version in use \($VERSION\) must\
    match the version found in the flagfile of the flowcell or run.\"
}

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    LIBRARY_PATH        \$root/lib
prepend-path    LD_LIBRARY_PATH     \$root/lib
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
