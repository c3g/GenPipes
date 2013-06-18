

###################
################### BLAST
###################
VERSION="2.2.28"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/blast/
mkdir -p $INSTALL_PATH
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$VERSION/ncbi-blast-$VERSION+-x64-linux.tar.gz"
tar -xvf ncbi-blast-$VERSION+-x64-linux.tar.gz
mv ncbi-blast-$VERSION+ $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC -  NCBI BLAST\"
}
module-whatis \" MUGQIC NCBI Blast: blast alignment\"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/blast/ncbi-blast-$VERSION+/bin
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/blast
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/blast/




