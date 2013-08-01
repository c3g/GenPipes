###################
################### MUGQIC TOOLS (hosted on svn for now)
###################
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools/tmp
cd $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools/tmp
VERSION="0.1"
git clone git@bitbucket.org:mugqic/mugqic_pipeline.git
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/mugqic_tools # where to install..
mkdir -p $INSTALL_PATH
cp -r mugqic_pipeline/tool_shed/* $INSTALL_PATH 
chmod -R 775 $INSTALL_PATH 

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - MUGQIC developped tools \"
}
module-whatis \"MUGQIC - MUGQIC developped tools \"
                       
set             root            \$::env(MUGQIC_INSTALL_HOME)/software/mugqic_tools
prepend-path    PATH            \$root/tools
prepend-path    PATH            \$root/perl-tools
prepend-path    PERL5LIB        \$root/perl-tools
setenv          R_TOOLS         \$root/R-tools
setenv          PERL_TOOLS      \$root/perl-tools

" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools
cd ..
rm -rf tmp
