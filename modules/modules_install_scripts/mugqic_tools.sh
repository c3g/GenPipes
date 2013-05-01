###################
################### MUGQIC TOOLS (hosted on svn for now)
###################
VERSION="0.1"
screen -S svn
ssh -l flefebvr  -L9443:esx-svn.genome.mcgill.ca:443 gallium.genome.mcgill.ca
# ctl+A +D, then enter when asking pwd
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/mugqic_tools # where to install..
mkdir -p $INSTALL_PATH
cp -r bioinformatics/R-tools bioinformatics/perl-tools bioinformatics/java-tools bioinformatics/tools $INSTALL_PATH 

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - MUGQIC developped tools \"
}
module-whatis \"MUGQIC - MUGQIC developped tools \"
                       
set             root            \$::env(MUGQIC_INSTALL_HOME)/software/mugqic_tools
prepend-path    PATH            \$root/tools
prepend-path    PATH            \$root/perl-tools
setenv          R_TOOLS         \$root/R-tools 

" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools


# mugqic/tools


