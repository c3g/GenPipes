###################
################### UCSC genome browser 'kent' bioinformatic utilities
###################
VERSION=`date +%Y%m%d`
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/ucsc/ucsc-$VERSION
INSTALL_DOWNLOAD=$MUGQIC_INSTALL_HOME/software/ucsc/tmp
mkdir -p $INSTALL_PATH $INSTALL_DOWNLOAD

# Download and extract
cd $INSTALL_DOWNLOAD
wget http://hgdownload.cse.ucsc.edu/admin/exe/userApps.src.tgz
tar zxvf userApps.src.tgz

# Compile
cd userApps
make
mv bin/* kentUtils.Documentation.txt $INSTALL_PATH
cd ..
chmod -R g+w $INSTALL_PATH


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tUCSC genome browser 'kent' bioinformatic utilities \"
}
module-whatis \"UCSC genome browser 'kent' bioinformatic utilities \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/ucsc/ucsc-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/ucsc/
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/ucsc/

rm -rf $INSTALL_DOWNLOAD
