

######### VirusFinder ##################
VERSION="031114"
NAME="art" # same could apply to all ucsc tools
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/art/
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download and extract 
wget http://www.niehs.nih.gov/research/resources/assets/docs/artsrcvanillaicecream${VERSION}linuxtgz.tgz
tar -xvzf artsrcvanillaicecream${VERSION}linuxtgz.tgz

# need module load gcc/4.7.2

#install gsl
cd art_src_VanillaIceCream_Linux/
wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
tar -xvzf gsl-1.16.tar.gz
mkdir $INSTALL_PATH/art_src_VanillaIceCream_Linux/gsl
cd gsl-1.16
./configure --prefix $INSTALL_PATH/art_src_VanillaIceCream_Linux/gsl
make
make install

export CFLAGS="$CFLAGS -I$INSTALL_PATH/art_src_VanillaIceCream_Linux/gsl/include" CPPFLAGS="$CPPFLAGS -I$INSTALL_PATH/art_src_VanillaIceCream_Linux/gsl/include" LDFLAGS="$LDFLAGS -L$INSTALL_PATH/art_src_VanillaIceCream_Linux/gsl/lib"

mkdir $INSTALL_PATH/art$VERSION
./configure --prefix $INSTALL_PATH/art$VERSION
make
make install

#TODO change /usr/bin/perl in /usr/bin/env perl in all perl script in $INSTALL_PATH/art$VERSION/bin
sed -i's/\/usr\/bin\/env/\/usr\/bin\/env perl/g' $INSTALL_PATH/art$VERSION/bin/aln2bed.pl
sed -i's/\/usr\/bin\/env/\/usr\/bin\/env perl/g' $INSTALL_PATH/art$VERSION/bin/combinedAvg.pl
sed -i's/\/usr\/bin\/env/\/usr\/bin\/env perl/g' $INSTALL_PATH/art$VERSION/bin/empDist.pl
sed -i's/\/usr\/bin\/env/\/usr\/bin\/env perl/g' $INSTALL_PATH/art$VERSION/bin/fastqReadAvg.pl
sed -i's/\/usr\/bin\/env/\/usr\/bin\/env perl/g' $INSTALL_PATH/art$VERSION/bin/map2bed.pl
sed -i's/\/usr\/bin\/env/\/usr\/bin\/env perl/g' $INSTALL_PATH/art$VERSION/bin/summation.pl

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - art \"
}
module-whatis \"MUGQIC - art \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME_DEV)/software/art/art$VERSION/bin/
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/art
mv .version $VERSION $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/art/



