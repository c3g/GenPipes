#!/bin/bash 
# Install software needed to run the ChIPSEQ pipeline
# Dependency of MUGQIC_HOME environment variable is assumed


main() {
  prologue
  InstallFreeBayes 
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

function isAvailable {
	###################
	################### isAvailable
	###################
	# Check if a module is available, matching by the module name, 
	# return the output of module avail command 
	moduleName=$1
	available=`module avail -l 2>&1 | grep $moduleName | sort -t "/" -k 3 -gr | sed -e 's/^.*\/+\('$moduleName'[/].*\) .$/\1/g' | head -1 | awk '{print $1}'`
	echo $available
}


function InstallFreeBayes()  {
	###################
	################### FreeBayes 
	###################
	# FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs, indels, MNPs (multi-nucleotide polymorphisms), 
	#and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.
	curDir=`pwd`;
	VERSION="v9.9.2-9-gfbf46fc"
	PACKAGE_NAME="FreeBayes"

	
	# 1 - Install software

	INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
	mkdir -p $INSTALL_PATH
	cd $INSTALL_PATH
	
	# A lot of complexity due to firewall restrictions on ABACUS :-p
	git clone https://github.com/ekg/freebayes.git
	cd freebayes
	https://github.com/pezmaster31/bamtools.github
	#Submodule 'bamtools' (git://github.com/pezmaster31/bamtools.git) registered for path 'bamtools'
	git clone https://github.com/pezmaster31/bamtools.git
  #Submodule 'intervaltree' (git://github.com/ekg/intervaltree.git) registered for path 'intervaltree'
  git clone https://github.com/ekg/intervaltree.git
  #Submodule 'vcflib' (git://github.com/ekg/vcflib.git) registered for path 'vcflib' and all its dependencies :-|
	git clone https://github.com/ekg/vcflib.git
	cd vcflib
	git clone https://github.com/ekg/fastahack.git
	git clone https://github.com/ekg/fsom.git
	git clone https://github.com/ekg/multichoose.git 
	git clone https://github.com/ekg/smithwaterman.git
	git clone https://github.com/ekg/tabixpp.git
	make CFLAGS="-O3 -D_FILE_OFFSET_BITS=64 -std=c++0x"
	make install CFLAGS="-O3 -D_FILE_OFFSET_BITS=64 -std=c++0x"
	cd $INSTALL_PATH/freebayes
  make CFLAGS="-O3 -D_FILE_OFFSET_BITS=64 -std=c++0x"
  make install CFLAGS="-O3 -D_FILE_OFFSET_BITS=64 -std=c++0x"
	
	# 2- Create module file
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 

	echo "#%Module1.0
	proc ModulesHelp { } {
	puts stderr \"\tMUGQIC - Adds Model-based Analysis of ChIP-Seq (MACS) tool to your environment \"
	}
	module-whatis \"MUGQIC - Adds Model-based Analysis of ChIP-Seq (MACS) tool to your environment \"
	set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"/freebayes
	prepend-path    PATH            \$root/bin
	setenv          FREEBAYES_BIN   \$root/bin 

	" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION

	cd $curDir
}

main
