
###################
################### R
###################

## Install R itself (libcairo must be installed?s)
VERSION="2.15.3"

# DEP_PATH is a URL or path to file with additional packages to be installed
DEP_PATH="http://bitbucket.org/mugqic/rpackages/raw/8e16c9322318a62ba74872504f0cef120803f1b7/DEPENDENCIES/modules_install_scripts/R_and_Bioconductor_packages.txt" 

# Download and compile and install
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/R/R-$VERSION # where to install..
mkdir -p $INSTALL_PATH
wget http://cran.r-project.org/src/base/R-${VERSION:0:1}/R-$VERSION.tar.gz
tar -xvf R-$VERSION.tar.gz
cd R-$VERSION
./configure --prefix=$INSTALL_PATH  # TEMP s--with-readline=yes --with-readline=no
make -j8
make install

## Install prefered add on packages (takes a loooong time)
wget $DEP_PATH
$INSTALL_PATH/bin/Rscript -e "source(\"http://bioconductor.org/biocLite.R\");deps=readLines(\"R_and_Bioconductor_packages.txt\");biocLite(deps,lib=.Library);biocLite(deps,lib=.Library)" # runs twice, sometimes mult attempts necessary

## chmod after installtion
chmod -R g+rw $INSTALL_PATH/lib64/R/library

# Module def file..
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Adds R to your environment \"
}
module-whatis \"MUGQIC - Adds R to your environment \"
                       
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/R/R-$VERSION
setenv          R_LIBS             \$root/lib64/R/library
prepend-path    MANPATH            \$root/share              
prepend-path    PATH               \$root/bin
prepend-path    LD_LIBRARY_PATH    \$root/lib64:/software/libraries/GotoBLAS_LAPACK/shared
#prepend-path   LD_LIBRARY_PATH    \$root/lib64:\$root/standalone:/software/libraries/GotoBLAS_LAPACK/shared
#prepend-path   CPATH              \$root/include

" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/R
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/R




