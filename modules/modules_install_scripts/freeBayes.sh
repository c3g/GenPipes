#!/bin/bash 
# Install software needed to run the Freebayes 
# Dependency of MUGQIC_HOME environment variable is assumed


main() {
  if [ ! -d "$MUGQIC_INSTALL_HOME" ];then
		prologue
	fi;
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
	PACKAGE_NAME="freeBayes"

	
	# 1 - Install software

	INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
	mkdir -p $INSTALL_PATH
	cd $INSTALL_PATH
	
	# change git config due to restrictions on ABACUS :-p
	git config --global url."https://".insteadOf git://
	git clone --recursive https://github.com/ekg/freebayes.git
	cd freebayes
	case $MUGQIC_INSTALL_HOME in
      	"/sb/programs/analyste")  
		make CFLAGS="-O3 -D_FILE_OFFSET_BITS=64 -std=c++0x"
		make install CFLAGS="-O3 -D_FILE_OFFSET_BITS=64 -std=c++0x"
	;;
	"/software/areas/genomics")
		module load cmake ifort_icc/12.0.4
		cat src/Makefile | sed -e 's/cmake/cmake -D CMAKE_C_COMPILER\=icc -D CMAKE\_CXX\_COMPILER\=icpc/g' | sed -e 's/CC\=g[+][+]/CC\=icpc/g' | sed -e 's/C\=gcc/\C\=icc/g'> src/tmpMake
 		mv src/tmpMake src/Makefile
		# some specific package's Makefiles need to be changed
		cat vcflib/tabixpp/Makefile | sed -e 's/CC\=gcc/CC\=icc/g' | sed -e 's/CPP\=g[+][+]/CPP\=icpc/g' > vcflib/tabixpp/tmpMake
		mv vcflib/tabixpp/tmpMake vcflib/tabixpp/Makefile		
		make CFLAGS="-O3 -D_FILE_OFFSET_BITS=64"
		cd vcflib/smithwaterman
		cat Makefile |  sed -e 's/CFLAGS\=/#CFLAGS\=/g' > tmpMake
		mv tmpMake Makefile
		make clean
		make CFLAGS="-O3 -D_FILE_OFFSET_BITS=64"
		cd $INSTALL_PATH/freebayes
		make
		make install
	;;
        *)  	
				make
				make install
    	;;
  	esac

	
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
