

###################
################### BLAST
###################
VERSION="2.2.28"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/blast/

# Set umask
umask 0002

mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$VERSION/ncbi-blast-$VERSION+-x64-linux.tar.gz"
tar -xvf ncbi-blast-$VERSION+-x64-linux.tar.gz
chmod -R g+w ncbi-blast-$VERSION+

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC -  NCBI BLAST\"
}
module-whatis \" MUGQIC NCBI Blast: blast alignment\"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/blast/ncbi-blast-$VERSION+/bin
prepend-path    PATH               \$root
" > $VERSION+

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION+\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/blast
mv .version $VERSION+ $MUGQIC_INSTALL_HOME/modulefiles/mugqic/blast/
rm ncbi-blast-$VERSION+-x64-linux.tar.gz

# Toxonomy DB
#  wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#  tar -xvf taxdb.tar.gz

