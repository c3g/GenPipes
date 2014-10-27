#!/bin/bash
# Install software needed to run the ChIPSEQ pipeline
# Dependency of MUGQIC_HOME environment variable is assumed


main() {
  #prologue
  #checkInstallWeblogo 
  checkInstallWeblogo282
  InstallHomer 
  InstallFindPeaks 
  InstallMACS
  InstallMACS2
}



function prologue {
	###################
	################### Prologue
	###################
	# Set the MUGQIC_INSTALL_HOME environment variable
	# Load MUGQIC module files
	umask 0002;
	## mugic module
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


function checkInstallNumpy {
	###################
	################### numpy 
	###################
	# Installing numpy on the loaded python version.
	# NumPy is an extension to the Python programming language, adding support for large, multi-dimensional arrays and matrices,
	# along with a large library of high-level mathematical functions to operate on these arrays. 
  
	python -c "import numpy" && numpyInstalled=1 || numpyInstalled=0;
	# If installed, keep the version of numpy, else install numpy 1.7.0
	if [ $numpyInstalled == 0 ]; then
		wget http://sourceforge.net/projects/numpy/files/NumPy/1.7.0/numpy-1.7.0.tar.gz/download -O numpy-1.7.0.tar.gz
		tar -xzvf numpy-1.7.0.tar.gz
		cd numpy-1.7.0
		python setup.py install
		cd ..
	fi;		
		# If numpy installation was no successful, a message will be send to the user
	python -c "import numpy";
}

function checkInstallGhostscript {
	###################
	################### Ghostscript
	################### Ghostscript is an interpreter for the PostScript language and for PDF.
  curDir=`pwd`;
  gs --version && gsInstalled=1 || gsInstalled=0;
  if [ $gsInstalled == 0 ]; then
			VERSION=9.02
			PACKAGE_NAME="ghostscript"
			INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
			mkdir -p $INSTALL_PATH
			mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/
			cd $INSTALL_PATH
			wget http://downloads.ghostscript.com/public/ghostscript-9.02.tar.gz -O ghostscript-9.02.tar.gz
			tar -xzvf ghostscript-9.02.tar.gz
			cd ghostscript-9.02
			./configure
			make
			make install
			
		  # 2- Create module file
	  	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 			
			echo "#%Module1.0
			proc ModulesHelp { } {
				puts stderr \"\tMUGQIC - Adds Ghostscript, an interpreter for the PostScript language and for PDF.\"
			}
			module-whatis \"MUGQIC - Adds Ghostscript, an interpreter for the PostScript language and for PDF. \"
 
			set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"
			prepend-path    PATH            \$root
			" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
	else
		  # 2- Create module file
      VERSION=`gs --version`
			PACKAGE_NAME="ghostscript"
			GHOSTSCRIPTPATH=`which gs | awk -F"/" 'BEGIN { OFS = "/"} {$NF=""; print $0}'`
	  	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 			
			echo "#%Module1.0
			proc ModulesHelp { } {
				puts stderr \"\tMUGQIC - Adds Ghostscript, an interpreter for the PostScript language and for PDF.\"
			}
			module-whatis \"MUGQIC - Adds Ghostscript, an interpreter for the PostScript language and for PDF. \"
 
			set             root            $GHOSTSCRIPTPATH
			prepend-path    PATH            \$root
			" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
	
	fi;
	
}

function checkInstallWeblogo {
	###################
	################### weblogo
	###################
	# Installing WebLogo
	# Dependencies
	# WebLogo version 3 is written in python. It is necessary to have python 2.5, 2.6 or 2.7 and the extension package numpy installed before WebLogo will run. 
	# WebLogo also requires a recent version of ghostscript to create PNG and PDF output, and pdf2svg to generate SVG output. 
		# 0 - Dependencies
	checkInstallNumpy; 
	checkInstallGhostscript; 
	moduleLoad=$( isAvailable "ghostscript" )
	if ["$moduleLoad"]; then
				module load $moduleLoad;
	fi;
	curDir=`pwd`
	VERSION="3.3"
	PACKAGE_NAME="weblogo"
		# 1 - Install software
	INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
	mkdir -p $INSTALL_PATH
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME

	cd $INSTALL_PATH

	module load mugqic/python/2.7.3
	wget http://weblogo.googlecode.com/files/weblogo-3.3.tar.gz -O weblogo-3.3.tar.gz
	tar -xzvf weblogo-3.3.tar.gz
	cd weblogo-3.3
	python setup.py install
	ln -s $INSTALL_PATH/weblogo $INSTALL_PATH/seqlogo

	echo "#%Module1.0
	proc ModulesHelp { } {
	puts stderr \"\tMUGQIC - Adds WebLogo, a tool for creating sequence logos from biological sequence alignments. \"
	}
	module-whatis \"MUGQIC -  Adds WebLogo, a tool for creating sequence logos from biological sequence alignments.  \"
	
	set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"
	setenv          WEBLOGO_HOME    \$root"/weblogo-3.3/"
	prepend-path    PATH            \$root"/weblogo-3.3/"
	" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
	# Default module version file
	echo "#%Module1.0
	set ModulesVersion \"$VERSION\"" > .version

	# Add permissions and install module
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME
	chmod -R ug+rwX $VERSION .version
	chmod -R o+rX $VERSION .version
	mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME

	# Clean up temporary installation files if any
	rm -rf $INSTALL_DOWNLOAD
	  cd $curDir
}

function checkInstallWeblogo282 {
	###################
	################### weblogo
	###################
	# Installing WebLogo v 2.8.2 from source files
	# Dependencies
	# WebLogo requires a recent version of ghostscript to create PNG and PDF output, and pdf2svg to generate SVG output. 
		# 0 - Dependencies
			checkInstallNumpy; 
	checkInstallGhostscript; 
	moduleLoad=$( isAvailable "ghostscript" )
	if ["$moduleLoad"]; then
				module load $moduleLoad;
	fi;
	curDir=`pwd`
	VERSION="2.8.2"
	PACKAGE_NAME="weblogo"
		# 1 - Install software
	INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
	mkdir -p $INSTALL_PATH
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME

	cd $INSTALL_PATH

	wget http://weblogo.berkeley.edu/release/weblogo.2.8.2.tar.gz -O weblogo.2.8.2.tar.gz
	tar -xzvf weblogo.2.8.2.tar.gz


	# Starting from guillimin Phase 2 (2014-02) you'll need to replace shebang #/usr/bin/perl with #/usr/bin/env perl and change the remaining parameters with the respective perl directives :-(
	# change /usr/bin/perl -w in /software/areas/genomics/phase2/software/weblogo/2.8.2/weblogo/seqlogo 
	# sed  -i "s/ perl -w/ perl \\nuse warnings;/" $MUGQIC_INSTALL_HOME/modulefiles/mugqic/weblogo/seqlogo 
	#	

 
	echo "#%Module1.0
	proc ModulesHelp { } {
	puts stderr \"\tMUGQIC - Adds WebLogo, a tool for creating sequence logos from biological sequence alignments. \"
	}
	module-whatis \"MUGQIC -  Adds WebLogo, a tool for creating sequence logos from biological sequence alignments.  \"
	
	set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"
	setenv          WEBLOGO_HOME    \$root/weblogo
	prepend-path    PATH            \$root/weblogo
	" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
	# Default module version file
	echo "#%Module1.0
	set ModulesVersion \"$VERSION\"" > .version

	# Add permissions and install module
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME
	chmod -R ug+rwX $VERSION .version
	chmod -R o+rX $VERSION .version
	mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME

	# Clean up temporary installation files if any
	rm -rf $INSTALL_DOWNLOAD
	cd $curDir
	
}

function InstallHomer {
	###################
	################### homer
	###################
	# Installing the basic HOMER software. HOMER will be installed in the same directory that you place the configureHomer.pl program.  
	# configureHomer.pl will attempt to check for required utilities and alert you to missing programs.
	# To install packages (Genomes), simply use the –install option and the name(s) of the package(s).
	# perl configureHomer.pl –install mouse (to download the mouse promoter set)
	# perl configureHomer.pl –install mm8    (to download the mm8 version of the mouse genome)
	# perl configureHomer.pl –install hg19r    (to download the hg19 repeat masked version of the human genome)
	curDir=`pwd`;
	VERSION="4.1"
	PACKAGE_NAME="homer"

	# 0 - Check and install dependencies
	module load $MUGQIC_INSTALL_HOME/modulefiles/mugqic/weblogo/3.3;
	moduleLoad=$( isAvailable "blat" )
	if [ "$moduleLoad" ]; then
		module load $moduleLoad;
	else
		echo "WARNING: Homer will be installed, but dependencies are not valid" $blat
	fi;

	# 1 - Install software
	INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
	mkdir -p $INSTALL_PATH
	cd $INSTALL_PATH
	wget http://biowhat.ucsd.edu/homer/configureHomer.pl 
	perl configureHomer.pl -install
	perl configureHomer.pl -install hg19
	perl configureHomer.pl -install mm10
	perl configureHomer.pl -install mm9
	perl configureHomer.pl -install rn5
	
	# Starting from guillimin Phase 2 (2014-02) you'll need to replace shebang #/usr/bin/perl with #/usr/bin/env perl and change the remaining parameters with the respective perl directives :-(
	#for fi in $(ls  $INSTALL_PATH/bin/*pl); 
	#do 
	#	cat $fi | sed -e 's/\/usr\/bin\/perl/\/usr\/bin\/env perl/g' > $fi.bak && mv $fi.bak $fi; 
	#	echo $fi
       	#        sed  -i "s/ perl -w/ perl \\nuse warnings;/" $fi
        #        sed -i "s/-I\/software\/areas\/genomics\/phase2\/software\/homer\/4.1\/.\/\/bin/\\nuse lib \'\/software\/areas\/genomics\/phase2\/software\/homer\/4.1\/bin\';/" $fi
	# done

	# 2- Create module file
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 

	echo "#%Module1.0
	proc ModulesHelp { } {
	puts stderr \"\tMUGQIC - Adds Homer, software for motif discovery and next generation sequencing analysis to your environment \"
	}
	module-whatis \"MUGQIC -  Adds Homer, software for motif discovery and next generation sequencing analysis to your environment \"
	set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"
	setenv          HOMER_HOME      \$root
	prepend-path    PATH            \$root/bin
	" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
	# Default module version file
	echo "#%Module1.0
	set ModulesVersion \"$VERSION\"" > .version
	
	# Add permissions and install module
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME
	chmod -R ug+rwX $VERSION .version
	chmod -R o+rX $VERSION .version
	mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME

	cd $curDir
}

function InstallFindPeaks()  {
	###################
	################### VancouverShortReads findPeaks 
	###################

	curDir=`pwd`
	VERSION="4.0.16"
	PACKAGE_NAME="VancouverShortReads"

	# 1 - Install software

	INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
	mkdir -p $INSTALL_PATH
	cd $INSTALL_PATH
	wget http://sourceforge.net/projects/vancouvershortr/files/latest/download -O VancouverShortR-4.0.16.tar.gz
	tar -xvzf VancouverShortR-4.0.16.tar.gz
	# jar file is now available on directory fp4

	# 2- Create module file
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 

	echo "#%Module1.0
	proc ModulesHelp { } {
	puts stderr \"\tMUGQIC - Adds Vancouver Short Reads FindPeaks tool to your environment \"
	}
	module-whatis \"MUGQIC - Adds Vancouver Short Reads FindPeaks tool to your environment \"
	
	set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"/fp4
	setenv          FINDPEAKS_JAR   \$root/FindPeaks.jar
	" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
	# Default module version file
	echo "#%Module1.0
	set ModulesVersion \"$VERSION\"" > .version
	
	# Add permissions and install module
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME
	chmod -R ug+rwX $VERSION .version
	chmod -R o+rX $VERSION .version
	mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME


	cd $curDir

}

function InstallMACS()  {
	###################
	################### MACS 
	###################
	# Model-based Analysis of ChIP-Seq (MACS), for identifying transcript factor binding sites.
	curDir=`pwd`;
	VERSION="1.4.2"
	PACKAGE_NAME="MACS"

	# This package depends on python 2.7.3
	module load mugqic/python/2.7.3

	# 1 - Install software

	INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
	mkdir -p $INSTALL_PATH
	cd $INSTALL_PATH
	wget --no-check-certificate https://github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz -O MACS-1.4.2-1.tar.gz
	tar -xvzf MACS-1.4.2-1.tar.gz

	cd MACS-1.4.2
	python setup.py install --prefix $INSTALL_PATH

	# 2- Create module file
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 

	echo "#%Module1.0
	proc ModulesHelp { } {
	puts stderr \"\tMUGQIC - Adds Model-based Analysis of ChIP-Seq (MACS) tool to your environment \"
	}
	module-whatis \"MUGQIC - Adds Model-based Analysis of ChIP-Seq (MACS) tool to your environment \"
	set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"
	prepend-path    PATH            \$root/bin
	prepend-path    PYTHONPATH      \$root/lib/python2.7/site-packages
	setenv          MACS_BIN        \$root/bin 
	setenv          MACS_LIB        \$root/lib/python2.7/site-packages

	" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
	# Default module version file
	echo "#%Module1.0
	set ModulesVersion \"$VERSION\"" > .version

	# Add permissions and install module
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME
	chmod -R ug+rwX $VERSION .version
	chmod -R o+rX $VERSION .version
	mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME

	cd $curDir
}

function InstallMACS2() {
	###################
	################### MACS2 
	###################
	# Model-based Analysis of ChIP-Seq (MACS), for identifying transcript factor binding sites.
	curDir=`pwd`;
	VERSION="2.0.10.09132012"
	PACKAGE_NAME="MACS"

	# This package depends on python 2.7.3
	# Prerequisites: Python version must be equal to *2.7* to run MACS. 
	# Numpy_ (>=1.3.0) are required to run MACS v2. 

	module load mugqic/python/2.7.3
	checkInstallNumpy
  
	# 1 - Install software
	INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
	mkdir -p $INSTALL_PATH
	cd $INSTALL_PATH
	wget --no-check-certificate https://pypi.python.org/packages/source/M/MACS2/MACS2-2.0.10.09132012.tar.gz -O MACS2-2.0.10.09132012.tar.gz
	tar -xvzf MACS2-2.0.10.09132012.tar.gz
	cd MACS2-2.0.10.09132012
	python setup.py install --prefix $INSTALL_PATH

	# 2- Create module file
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 

	echo "#%Module1.0
	proc ModulesHelp { } {
	puts stderr \"\tMUGQIC - Adds Model-based Analysis of ChIP-Seq (MACS) tool to your environment \"
	}
	module-whatis \"MUGQIC - Adds Model-based Analysis of ChIP-Seq (MACS) tool to your environment \"
	set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"
	prepend-path    PATH            \$root/bin
	prepend-path    PYTHONPATH      \$root/lib/python2.7/site-packages
	setenv          MACS_BIN        \$root/bin 
	setenv          MACS_LIB        \$root/lib/python2.7/site-packages

	" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
	# Default module version file
	echo "#%Module1.0
	set ModulesVersion \"$VERSION\"" > .version
	
	# Add permissions and install module
	mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME
	chmod -R ug+rwX $VERSION .version
	chmod -R o+rX $VERSION .version
	mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME

	cd $curDir
}


main
