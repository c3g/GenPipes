# Install software needed to run the ChIPSEQ pipeline
# Dependency of MUGQIC_HOME.. is assumed

#### Epilogue
umask 0002
## mugic module
HOST=`hostname`;
DNSDOMAIN=`dnsdomainname`;
if [[ $HOST == abacus* ]]; then
 export MUGQIC_INSTALL_HOME=/sb/programs/analyste
elif [[ $HOST == lg-* ]]; then
 export MUGQIC_INSTALL_HOME=/software/areas/genomics
elif [[ $HOST == ip03 ]]; then
 export MUGQIC_INSTALL_HOME=/mnt/scratch_mp2/bourque/bourque_group/opt
elif [[ $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then
 export MUGQIC_INSTALL_HOME=/sb/programs/analyste
elif [[ $DNSDOMAIN == guillimin.clumeq.ca ]]; then
 export MUGQIC_INSTALL_HOME=/software/areas/genomics
elif [[ $DNSDOMAIN == m ]]; then
 export MUGQIC_INSTALL_HOME=/mnt/scratch_mp2/bourque/bourque_group/opt
fi
module use $MUGQIC_INSTALL_HOME/modulefiles



###################
################### homer
###################
# Installing the basic HOMER software. HOMER will be installed in the same directory that you place the configureHomer.pl program.  
# configureHomer.pl will attempt to check for required utilities and alert you to missing programs.
# To install packages (Genomes), simply use the –install option and the name(s) of the package(s).
# perl configureHomer.pl –install mouse (to download the mouse promoter set)
# perl configureHomer.pl –install mm8    (to download the mm8 version of the mouse genome)
# perl configureHomer.pl –install hg19r    (to download the hg19 repeat masked version of the human genome)

curDir=`pwd`
VERSION="4.1"
PACKAGE_NAME="homer"

# 1 - Install software
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
wget http://biowhat.ucsd.edu/homer/configureHomer.pl 
perl configureHomer.pl -install

# 2- Create module file
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/ 

echo "#%Module1.0
proc ModulesHelp { } {
 puts stderr \"\tMUGQIC - Adds Homer, software for motif discovery and next generation sequencing analysis to your environment \"
}
module-whatis \"MUGQIC -  Adds Homer, software for motif discovery and next generation sequencing analysis to your environment \"
 
set             root            \$::env(MUGQIC_INSTALL_HOME)/software/"$PACKAGE_NAME"/"$VERSION"
setenv          HOMER_HOME      \$root/bin
prepend-path    PATH            \$root/bin
" > $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$PACKAGE_NAME/$VERSION
cd $curDir


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
wget http://sourceforge.net/projects/vancouvershortr/files/latest/download
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

cd $curDir


###################
################### MACS 
###################

curDir=`pwd`
VERSION="1.4.2"
PACKAGE_NAME="MACS"

# 1 - Install software

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$PACKAGE_NAME/$VERSION # where to install..
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
wget https://github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz
tar -xvzf MACS-1.4.2-1.tar.gz

# This package depends on python 2.7.3
module load mugqic/python/2.7.3
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

cd $curDir
