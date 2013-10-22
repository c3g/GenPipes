###################
################### MUGQIC TOOLS (hosted on svn for now)
###################
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools/tmp $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools/tmp/untar/
cd $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools/tmp
VERSION="1.0"
wget https://bitbucket.org/mugqic/mugqic_tools/get/v${VERSION}.tar.gz
tar -xvf v${VERSION}.tar.gz -C untar/
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/mugqic_tools # where to install..
ARCHIVE_PATH=$MUGQIC_INSTALL_HOME/archive/mugqic_tools
mkdir -p $INSTALL_PATH
cp untar/mugqic-mugqic_tools-*/*  $INSTALL_PATH
chmod -R 775 $INSTALL_PATH 
mv v${VERSION}.tar.gz $ARCHIVE_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - MUGQIC developped tools \"
}
module-whatis \"MUGQIC - MUGQIC developped tools \"
                       
set             root            \$::env(MUGQIC_INSTALL_HOME)/software/mugqic_tools
prepend-path    PATH            \$root/tools
prepend-path    PATH            \$root/perl-tools
prepend-path    PATH            \$root/python-tools
prepend-path    PERL5LIB        \$root/perl-tools
setenv          R_TOOLS         \$root/R-tools
setenv          PERL_TOOLS      \$root/perl-tools
setenv          PYTHON_TOOLS    \$root/python-tools

" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools
cd ..
rm -rf tmp
