#!/bin/bash


### NOTE (FL): rnmmer is a pain to install:
# - requires 2.x version of hmmer
# - requires XML::Simple module
# - see http://blog.karinlag.no/2013/10/rnammer-install/ for detauls


# Installation notes:
# RNAMMER requires the older version of hmmsearch (v2).
# You can obtain the hmmsearch_v2 at ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz[ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz].
# After building the software, we suggest you rename this version of hmmsearch as 'hmmsearch2'.
# In the 'rnammer' software configuration, edit the rnammer script to point
# $HMMSEARCH_BINARY = "/path/to/hmmsearch2";
# Be sure that rnammer functions correctly by executing it on their provided sample data.  RNAMMER is quite useful, but the current implementation is not robust to error, so check carefully.
# 
# 


# RNAmmer 1.2		INSTALLATION INSTRUCTIONS
# 
# 
# DESCRIPTION
# 
# RNAmmer 1.2 predicts 5s/8s, 16s/18s, and 23s/28s ribosomal RNA  in full genome
# sequences. The method is described in detail in the following article: 
# 
# RNammer: consistent annotation of rRNA genes in genomic sequences.
# Lagesen K, Hallin PF, Roedland E, Staerfeldt HH, Rognes T Ussery DW.
# Nucleic Acids Res. Apr 22, 2007.
# 
# More information about the method can be found at:
# 
# 	http://www.cbs.dtu.dk/services/RNAmmer/
# 	http://www.cbs.dtu.dk/ws/RNAmmer/
# 
# DOWNLOAD
# 
# The RNAmmer package  is a property of Center for Biological Sequence Analysis.
# It may be downloaded only by special agreement.  For academic users there is a
# download site at:
# 
# 	http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?rnammer
# 
# Other users are requested to contact software@cbs.dtu.dk.
# 
# PRE-INSTALLATION
# 
# RNAmmer will run on the most common UNIX platforms e.g.  Linux, SunOS etc.
# Make sure that your server has a complete UNIX installation. Specifically, the
# following programs are necessary:
# 
# 	Perl scripting language with support for Getopt::Long module installed
# 	hmmsearch	- the profile HMM search program from HMMER
# 
# INSTALLATION
# 
# 1. Decide where you wish to keep the software. Uncompress and untar
#    the package in that location:
# 
# 	cat rnammer-1.2.tar.Z | gunzip | tar xvf -
# 
#    This will produce a directory  'rnammer-1.2.'.
# 
# 2. Edit path specifications in the program file 'rnammer'.
# 
# 3. Test RNAmmer on the test sequences shipped with the package:
# 
# 	perl rnammer -S bac -m lsu,ssu,tsu -gff - < example/ecoli.fsa
# 
# 4. Move or copy the 'rnammer' script to a directory in the users' path.
# 
# 5. Move or copy the 'man/rnammer.1' file to a appropriate location  in your manual
#    system. 
# 
# 6. Enjoy ...
# 
# 
# PROBLEMS AND QUESTIONS
# 
# In case of technical problems (bugs etc.) please contact packages@cbs.dtu.dk.
# 
# Questions on the scientific aspects of the RNAmmer method  should go to 
# Peter Hallin pfh@cbs.dtu.dk.
# 
# 
# CBS, July 19 2007





# rm -rf /mnt/parallel_scratch_mp2_wipe_on_august_2014/bourque/bourque_group/analyste/software/rnammers

# HMMSEARCH2PATH
module load mugqic/hmmer/2.3.2
HMMSEARCH2PATH=`which hmmsearch`

# Required perl installation
module load mugqic/perl/5.18.2
# cpan install XML::Simple

# need to ln -s  hmmmer2
SOFTWARE=rnammer
VERSION=1.2
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE/$SOFTWARE-$VERSION 
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
wget https://www.dropbox.com/s/5j05i5vs2s18tee/rnammer-$VERSION.src.tar.Z -O rnammer-$VERSION.src.tar.Z
tar xvfz rnammer-$VERSION.src.tar.Z
mv -if rnammer-$VERSION.src.tar.Z $MUGQIC_INSTALL_HOME/archive


# link to hmmsearch 2.3.2
ln -s $HMMSEARCH2PATH hmmsearch2
sed -i "s,\$HMMSEARCH_BINARY = \"/usr/cbs/bio/bin/linux64/hmmsearch\";,\$HMMSEARCH_BINARY = \"$PWD/hmmsearch2\";,g" rnammer 

# Paths
sed -i "s,\"/usr/cbs/bio/src/rnammer-1.2\";,\"$PWD\";,g" rnammer 

# perl
sed -i s,"#\!/usr/bin/perl,#\!/usr/bin/env perl,g" rnammer 
sed -i s,"#\!/usr/bin/perl,#\!/usr/bin/env perl,g" xml2fsa
sed -i s,"#\!/usr/bin/perl,#\!/usr/bin/env perl,g" xml2gff
sed -i s,"#\!/usr/bin/perl,#\!/usr/bin/env perl,g" core-rnammer

# Crap -cpu 1 issue
# http://blog.karinlag.no/2013/10/rnammer-install/
sed -i "s,--cpu 1,,g" core-rnammer

## Test
perl rnammer -S bac -m lsu,ssu,tsu -gff - < example/ecoli.fsa

# Add permissions and install software
cd ..
chmod -R ug+rwX .
chmod -R o+rX .


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ; 
}
module-whatis \"$SOFTWARE  \" ; 

set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;  
prepend-path    PATH                \$root ;
" > $VERSION


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


