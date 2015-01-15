#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

### NOTE (FL): rnmmer is a pain to install:
# - requires 2.x version of hmmer
# - requires XML::Simple module
# - see http://blog.karinlag.no/2013/10/rnammer-install/ for details


# Installation notes:
# RNAMMER requires the older version of hmmsearch (v2).
# You can obtain the hmmsearch_v2 at ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz[ftp://selab.janelia.org/pub/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz].
# After building the software, we suggest you rename this version of hmmsearch as 'hmmsearch2'.
# In the 'rnammer' software configuration, edit the rnammer script to point
# $HMMSEARCH_BINARY = "/path/to/hmmsearch2";
# Be sure that rnammer functions correctly by executing it on their provided sample data.  RNAMMER is quite useful, but the current implementation is not robust to error, so check carefully.


# RNAmmer 1.2   INSTALLATION INSTRUCTIONS
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
#   http://www.cbs.dtu.dk/services/RNAmmer/
#   http://www.cbs.dtu.dk/ws/RNAmmer/
# 
# DOWNLOAD
# 
# The RNAmmer package  is a property of Center for Biological Sequence Analysis.
# It may be downloaded only by special agreement.  For academic users there is a
# download site at:
# 
#   http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?rnammer
# 
# Other users are requested to contact software@cbs.dtu.dk.
# 
# PRE-INSTALLATION
# 
# RNAmmer will run on the most common UNIX platforms e.g.  Linux, SunOS etc.
# Make sure that your server has a complete UNIX installation. Specifically, the
# following programs are necessary:
# 
#   Perl scripting language with support for Getopt::Long module installed
#   hmmsearch - the profile HMM search program from HMMER
# 
# INSTALLATION
# 
# 1. Decide where you wish to keep the software. Uncompress and untar
#    the package in that location:
# 
#   cat rnammer-1.2.tar.Z | gunzip | tar xvf -
# 
#    This will produce a directory  'rnammer-1.2.'.
# 
# 2. Edit path specifications in the program file 'rnammer'.
# 
# 3. Test RNAmmer on the test sequences shipped with the package:
# 
#   perl rnammer -S bac -m lsu,ssu,tsu -gff - < example/ecoli.fsa
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

SOFTWARE=rnammer
VERSION=1.2
ARCHIVE=$SOFTWARE-$VERSION.src.tar.Z
# Must send a request for academic users to get this URL
ARCHIVE_URL=https://www.dropbox.com/s/5j05i5vs2s18tee/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

# Specific commands to extract and build the software
# $INSTALL_DIR and $INSTALL_DOWNLOAD have been set automatically
# $ARCHIVE has been downloaded in $INSTALL_DOWNLOAD
build() {
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  # Archive has no root dir; therefore extract it directly in destination dir
  tar zxvf $INSTALL_DOWNLOAD/$ARCHIVE --directory=$INSTALL_DIR/$SOFTWARE_DIR

  # Various custom modifications due to poor original coding...
  cd $INSTALL_DIR/$SOFTWARE_DIR
  # Remove hmmsearch absolute path (IMPORTANT: module mugqic/hmmer/2.3.2 must be loaded to run rnammer!)
  sed -i "s,HMMSEARCH_BINARY = \"\/[^ ]*hmmsearch,HMMSEARCH_BINARY = \"hmmsearch,g" rnammer
  # Also replace hmmsearch file check with 'which' to avoid error when using hmmsearch via module mugqic/hmmer/2.3.2
  sed -i "s,err_exit ( \"read_configuration(): Binary 'hmmsearch' not found (specify with '\[hmmsearch\]=)\" ) unless -f \$config{hmmsearch},use File::Which;\n\terr_exit ( \"read_configuration(): Binary 'hmmsearch' not found (specify with '[hmmsearch]=)\" ) unless which(\$config{hmmsearch}),g" core-rnammer

  # Remove perl absolute path (IMPORTANT: module mugqic/perl with XML::Simple must be loaded to run rnammer!)
  sed -i "s,PERL = \"\/[^ ]*perl,PERL = \"perl,g" rnammer

  # Update rnammer program path with dynamic path cmd
  sed -i "s,my \$INSTALL_PATH = \"/usr/cbs/bio/src/rnammer-1.2\",use FindBin;\nmy \$INSTALL_PATH = \$FindBin::Bin," rnammer

  # Update Perl script shebangs
  sed -i s,"#\!/usr/bin/perl,#\!/usr/bin/env perl,g" core-rnammer rnammer xml2fsa xml2gff

  # Crap -cpu 1 issue
  # Cf: http://blog.karinlag.no/2013/10/rnammer-install/
  sed -i "s,--cpu 1,,g" core-rnammer

  # Test RNAmmer
  set +e
  echo "Testing RNAmmer..."
  module load mugqic/hmmer/2.3.2
  module load mugqic/perl/5.18.2
  ./rnammer -S bac -m lsu,ssu,tsu -xml ecoli.xml -gff ecoli.gff -h ecoli.hmmreport < example/ecoli.fsa
  echo "The output files ecoli.xml, ecoli.gff and ecoli.hmmreport should be identical to the corresponding files in '$SOFTWARE_DIR/example'"
  echo "diff example/ecoli.xml ecoli.xml"
  diff example/ecoli.xml ecoli.xml
  echo "diff example/ecoli.gff ecoli.gff"
  diff example/ecoli.gff ecoli.gff
  echo "diff example/ecoli.hmmreport ecoli.hmmreport"
  diff example/ecoli.hmmreport ecoli.hmmreport
  echo "Testing RNAmmer finished."
  set -e
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

prereq                              mugqic/hmmer/2.3.2
set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
