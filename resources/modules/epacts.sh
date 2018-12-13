#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

# NOTES:
# - The script 	$INSTALL_DIR/$SOFTWARE_DIR/bin/epacts download downloads ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz with perl NET::FTP, this fails on abacus. If we ever want to patch:

# Bon, oui il se trouve que ftp soit problèmatique sur abacus car il utilise des ports tcp aléatoires.  Heureusement, c'est assez peu utilisé normalement.  La pluspart des téléchargement utilisent http, https, sftp ou encore des protocole spécialisés comme UDT ou EGA.
#
# Alors j'ai solutionné le problème comme suit:
#
# 1- J'ai constaté que le serveur distribue les mêmes données sur http aussi (pas toujours le cas mais souvent):
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/
#
# 2- J'ai patché le programme afin qu'il récupère les données sur http plus tôt que ftp.  Je l'ai testé, les fichiers sont sur abacus maintenant.
#
# La nouvelle commande ainsi créée, epacts-download-http, peut télécharger les fichiers en question.  Bien sûr, j'ai laissé l'originale intacte
#
# Marc-andré
#
# P.s: Voici le diff qui montre le changement que j'ai fait.
#
# diff -u epacts-download epacts-download-http
# --- epacts-download     2015-07-29 14:50:24.770927000 -0400
# +++ epacts-download-http        2015-07-29 16:20:54.238175000 -0400
# @@ -6,6 +6,7 @@
#  use FindBin;
#  use lib "$FindBin::Bin";
#  use Net::FTP;
# +use LWP::Simple;
#
#  my $man = 0;
#  my $help = 0;
# @@ -26,7 +27,11 @@
#  my $dir = "/vol1/ftp/technical/reference";
#  my $fasta = "human_g1k_v37.fasta.gz";
#  my $fai = "human_g1k_v37.fasta.fai";
# +my $result = 500;
# +my $url = '';
#
# +# ftp won't work falling back to http
# +=pod
#  print "Connecting to $hostname\n";
#  my $ftp = Net::FTP->new("$hostname", Debug => 0) or die "Cannot connect to $hostname $@";
#  $ftp->login("anonymous",'-anonymous@') || die "Cannot login ", $ftp->message;
# @@ -47,6 +52,26 @@
#  print "Moving $fasta to $datadir/\n";
#  rename("$fasta","$datadir/$fasta");
#  $ftp->quit;
# +=cut
# +
# +print "Downloading $fai..\n";
# +$url = "http://$hostname$dir/$fai";
# +$result = getstore($url, $fai);
# +if ( $result != 200 ) {
# +       die "Download of $fai failed, url: $url, result: $result";
# +}
# +print "Moving $fai to $datadir/\n";
# +rename("$fai","$datadir/$fai");
# +
# +print "Downloading $fasta..\n";
# +$url = "http://$hostname$dir/$fasta";
# +getstore($url,$fasta);
# +if ( $result != 200 ) {
# +       die "Download of $fasta failed, url: $url, result: $result";
# +}
# +print "Moving $fasta to $datadir/\n";
# +rename("$fasta","$datadir/$fasta");
# +
#
#  print "Decompressing the file\n";
#  my $cmd = "gzip -d $datadir/$fasta";


# - http://csg.sph.umich.edu/kang/epacts/download/$ARCHIVE was originally blocketd on abacus, now fixed.


# (See http://genome.sph.umich.edu/wiki/EPACTS for comprehensive documentation)
#
# == EPACTS Installation Details  ==
#
# If you want to use EPACTS in an Ubuntu platform, following the step below
#
# * Download EPACTS source distribution at http://www.sph.umich.edu/csg/kang/epacts/download/EPACTS-3.0.0.tar.gz (99MB)
# * Uncompress EPACTS package, and install the package using the following set of commands
#
#   tar xzvf EPACTS-$(VERSION).tar.gz
#   cd EPACTS-$(VERSION)
#   ./configure --prefix=/path/to/install
#   (If you have libraries in non-standard directory, please try
#    ./configure --prefix=/path/to/install LDFLAGS=-L/path/to/library to include the directory containing libR.so)
#   make -j12
#   make install
#
# (Important Note: '''make sure to specify --prefix=/path/to/install''' to avoid installing to the default path /usr/local/, which you may not have the permission. /home/your_userid/epacts might be a good one, if you are not sure where to install)
#
# * Now ${EPACTS_DIR} represents the '/path/to/install' directory
#
# * Download the reference FASTA files from 1000 Genomes FTP automatically by running the following commands
#
#   ${EPACTS_DIR}/bin/epacts download
#
#  (For advanced users, to save time for downloading the FASTA files (~900MB), you may copy a local copy of GRCh37 FASTA file and the index file to ${EPACTS_DIR}/share/EPACTS/)
#
# *Perform a test run by running the following command
#
#   ${EPACTS_DIR}/bin/test_run_epacts.sh


PREREQS="mugqic/gnuplot/4.6.6 mugqic/R_Bioconductor/3.1.2_3.0"
SOFTWARE="EPACTS"
VERSION="3.2.6"  
ARCHIVE="$SOFTWARE-$VERSION.tar.gz"  
ARCHIVE_URL="http://csg.sph.umich.edu/kang/epacts/download/$ARCHIVE"
SOFTWARE_DIR=$SOFTWARE-$VERSION 

build() {
  cd $INSTALL_DOWNLOAD
  tar zxvf $ARCHIVE 
  cd $SOFTWARE_DIR
	module load $PREREQS
  ./configure --prefix=$INSTALL_DIR/$SOFTWARE_DIR  
  make -j12 
  make install  
	
	$INSTALL_DIR/$SOFTWARE_DIR/bin/epacts download
	$INSTALL_DIR/$SOFTWARE_DIR/bin/test_run_epacts.sh &> $INSTALL_DIR/$SOFTWARE_DIR/test_run_epacts.sh.out   # $HOME/test/bin/test_run_epacts.sh
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

prereq mugqic/gnuplot/4.6.6 mugqic/R_Bioconductor/3.1.2_3.0

set             root                $INSTALL_DIR/$SOFTWARE_DIR ;
prepend-path    PATH                \$root/bin ;  
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
