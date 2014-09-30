
###################
################### STAR
###################
VERSION="2.4.0d"
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/star
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download and extract
wget http://github.com/alexdobin/STAR/archive/STAR_2.4.0d.tar.gz -O STAR_2.4.0d.tar.gz
tar -xvf STAR_$VERSION.tar.gz
rm  STAR_$VERSION.tar.gz

##INSTALL
cd STAR-STAR_$VERSION
make -j8
cd ..
chmod -R ug+rwX STAR-STAR_$VERSION
chmod -R o+rX STAR-STAR_$VERSION


# Module file
echo "#%Module1.0
proc ModulesHelp { } {STAR RNA-aligner \"
}
module-whatis \"MUGQIC - star \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME_DEV)/software/star/STAR-STAR_$VERSION
prepend-path    PATH               \$root
" > $VERSION



# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/star
mv .version $VERSION $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/star/

