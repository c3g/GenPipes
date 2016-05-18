#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

################################################################################
# This is a module install script template which should be copied and used for
# consistency between module paths, permissions, etc.
# Only lines marked as "## TO BE ADDED/MODIFIED" should be, indeed, modified.
# Also, once modified, delete this commented-out header and the ## comments
################################################################################

HOST=`hostname`;

DNSDOMAIN=`dnsdomainname`;

SOFTWARE=VarDictJava 
VERSION=1.4.5
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
ARCHIVE_URL=https://github.com/AstraZeneca-NGS/VarDictJava/archive/v$VERSION.tar.gz  
SOFTWARE_DIR=$SOFTWARE-$VERSION 

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE  

  ##Unzip compile version
  cd $SOFTWARE_DIR/dist
  unzip VarDict-${VERSION}.zip
  mv VarDict-${VERSION}/* ../
  
  ##Modify Vardict bash script
  cd ../
  if [[ $HOST == abacus* || $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then

     sed -i.bak -e "s/^DEFAULT_JVM_OPTS='/DEFAULT_JVM_OPTS=\'\"-Djava.io.tmpdir=\/lb\/scratch/\" \"-XX:ParallelGCThreads=1\" \"-Dsamjdk.buffer_size=4194304\" /g" bin/VarDict

  elif [[ $HOST == lg-* || $DNSDOMAIN == guillimin.clumeq.ca ]]; then

     sed -i.bak -e "s/^DEFAULT_JVM_OPTS='/DEFAULT_JVM_OPTS=\'\"-Djava.io.tmpdir=\/localscratch\/\" \"-XX:ParallelGCThreads=1\" \"-Dsamjdk.buffer_size=4194304\" /g" bin/VarDict

  elif [[ $BQMAMMOUTH == "mp2" ]]; then

     sed -i.bak -e "s/^DEFAULT_JVM_OPTS='/DEFAULT_JVM_OPTS=\'\"-Djava.io.tmpdir=\/lb\/scratch\/\" \"-XX:ParallelGCThreads=1\" \"-Dsamjdk.buffer_size=4194304\" /g" bin/VarDict

  fi  
  
  ##Install Vardict perl support files  
  rm -R VarDict
  git clone https://github.com/AstraZeneca-NGS/VarDict.git
 
  # Install software
  cd $INSTALL_DOWNLOAD
  mv -i $SOFTWARE_DIR $INSTALL_DIR/
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin ; 
setenv          VARDICT_HOME        \$root ;
setenv          VARDICT_BIN         \$root/VarDict ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
