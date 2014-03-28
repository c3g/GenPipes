#!/bin/bash

# NOTE: 
# - Perl module DB_File is a dependency for the Transdecoder part of Trinotate. This module depends on some version BerkeleyDB which was not present on MAmmouth...
# - 


#TEMPDIR=`mktemp -d -t $me.XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX` && cd $TEMPDIR 
#echo "Working in $TEMPDIR"
# rm -rf /mnt/parallel_scratch_mp2_wipe_on_august_2014/bourque/bourque_group/analyste/software/trinotate
SOFTWARE=trinotate  
VERSION=20131110 
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
wget http://downloads.sourceforge.net/project/trinotate/Trinotate_r$VERSION.tar.gz
tar xvf Trinotate_r$VERSION.tar.gz 
#cd Trinotate_r$VERSION 
wget "http://sourceforge.net/projects/trinotate/files/TRINOTATE_RESOURCES/TrinotateSqlite.sprot.$VERSION.db.gz/download" -O Trinotate_r$VERSION/Trinotate.sqlite.gz
gunzip Trinotate_r$VERSION/Trinotate.sqlite.gz

## Uniprot, PFAM : see install scripts for those


# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .
mv -if Trinotate_r$VERSION.tar.gz $MUGQIC_INSTALL_HOME/archive  

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ; 
}
module-whatis \"$SOFTWARE  \" ; 
prereq                               mugqic/trinity/$VERSION
#prereq                              mugqic/sqlite3
#prereq                               mugqic/blast/2.2.29+
#prereq                               mugqic/hmmer/3.1b1
#prereq                               mugqic/signalp/4.1
#prereq                               mugqic/tmhmm/2.0c
#prereq                               mugqic/rnammer/1.2
#prereq                               mugqic/cd-hit/4.5.4-2011-03-07
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/Trinotate_r$VERSION ;
prepend-path    PATH                \$root ; 
setenv          TRINOTATE_HOME         \$root;
setenv         TRINOTATE_SQLITE       \$root/Trinotate.sqlite;
" > $VERSION

# NOTES:
# - Assuming trinotate and trinity version follow one another
# - Assuming sqlite is already available on the system

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE

