#!/bin/bash
# Install software needed to run the Freebayes 
# Dependency of MUGQIC_HOME environment variable is assumed


main() {
  if [ ! -d "$MUGQIC_INSTALL_HOME" ];then
                prologue
        fi;
  Install
}


function prologue {
  # This script runs operations common to all installtions, such as setting the mask and defining $MUGQIC_INSTALL_HOME. It should also be added to the user's bash profile

  # Set umask
  umask 0002

  # Guess cluster and set $MUGQIC_INSTALL_HOME accordingly
  HOST=`hostname`;
  DNSDOMAIN=`dnsdomainname`;
  if [[ $HOST == abacus* ]]; then
   export MUGQIC_INSTALL_HOME=/sb/programs/analyste
  elif [[ $HOST == lg-* ]]; then
   export MUGQIC_INSTALL_HOME=/software/areas/genomics
  elif [[ $HOST == ip03 ]]; then
   export MUGQIC_INSTALL_HOME=/mnt/lustre03/bourque/bourque_group/opt
  elif [[ $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then
   export MUGQIC_INSTALL_HOME=/sb/programs/analyste
  elif [[ $DNSDOMAIN == guillimin.clumeq.ca ]]; then
   export MUGQIC_INSTALL_HOME=/software/areas/genomics
  elif [[ $DNSDOMAIN == m ]]; then
     export MUGQIC_INSTALL_HOME=/mnt/lustre03/bourque/bourque_group/opt
  fi

  # Module use
  module use $MUGQIC_INSTALL_HOME/modulefiles
}

function Install()  {
        ###################
        ################### The Gem library
        ###################
        # The GEM mapper, The GEM RNA mapper, The GEM mappability , and others
        curDir=`pwd`;
        VERSION="v1.315"
        PACKAGE_NAME="gemLibrary"

        
        # 1- Install software

        INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
        mkdir -p $INSTALL_PATH
        cd $INSTALL_PATH
        wget http://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download
        wget http://downloads.sourceforge.net/project/gemlibrary/gem-library/Binary%20pre-release%202/GEM-binaries-Linux-x86_64-core_2-20121106-022124.tbz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fgemlibrary%2Ffiles%2Fgem-library%2FBinary%2520pre-release%25202%2F&ts=1377722916&use_mirror=softlayer-dal
        tar -xvf GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2
        tar -xvf GEM-binaries-Linux-x86_64-core_2-20121106-022124.tbz2 
        chmod -R g+w GEM-binaries-Linux-x86_64-core_2-20121106-022124 GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin
        rm -rf GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2 GEM-binaries-Linux-x86_64-core_2-20121106-022124.tbz2

        # 2- Create module file
        mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 

        echo "#%Module1.0
        proc ModulesHelp { } {
        puts stderr \"\tMUGQIC - Adds  The Gem library (gem) tool to your environment \"
        }
        module-whatis \"MUGQIC - Adds the Gem library (The GEM mapper, The GEM RNA mapper, The GEM mappability , and others) tool to your environment \"
        set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"
        prepend-path    PATH            \$root/GEM-binaries-Linux-x86_64-core_2-20121106-022124
        prepend-path    PATH            \$root/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin

        " > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION

        # 3- Version file
        echo "#%Module1.0
        set ModulesVersion \"$VERSION\"
        " > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/.version

        cd $curDir
}

main
