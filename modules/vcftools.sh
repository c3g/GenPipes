
###################
################### Vcftools
###################
VERSION="0.1.11"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/vcftools/
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_PATH/archive $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, compile
wget http://sourceforge.net/projects/vcftools/files/vcftools_$VERSION.tar.gz/download
tar zxvf vcftools_$VERSION.tar.gz
cd vcftools_$VERSION
make
cd ..

# Move to install path
mv vcftools_$VERSION $INSTALL_PATH
chmod -R g+w $INSTALL_PATH/vcftools_$VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - vcftools \"
}
module-whatis \"MUGQIC - vcftools \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/vcftools/vcftools_$VERSION
prepend-path    PATH               \$root/bin
prepend-path    PERL5LIB           \$root/lib/perl5/site_perl
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/vcftools
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/vcftools
mv $INSTALL_DOWNLOAD/vcftools_$VERSION.tar.gz $INSTALL_PATH/archive/
rm -rf $INSTALL_DOWNLOAD
