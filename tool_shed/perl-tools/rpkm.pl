#!/usr/bin/perl
use Time::HiRes qw( time );

use strict;
use warnings;
use Getopt::Std;

use vars qw($VERSION);
$VERSION = '1.0';

&main;

sub main {

  my %opts = ( r => undef );
  getopts( 'r:', \%opts );

  if ( !defined( $opts{'r'} ) ) {

    die(
    qq/
Usage:   rpkm.pl [options] -r <Velvet\/Oases run directory> -o <output>
Version: $VERSION
Options: -r        Velvet\/Oases run directory (Directory with LastGraph, Sequences, contig-ordering, Log files).
\n/
    )
  }

  my $readsCountFileName = $opts{'r'};

  open( READS_FILE, "$readsCountFileName" ) or die "Can't open $readsCountFileName\n";
  $_ = <READS_FILE>;
  chomp;
  my (undef, @samplesTotalNbReads) = split(/\t/);
  my @nbReadsPerMillion;
  for my $totalNbReads (@samplesTotalNbReads) {
    my $readsPerMillion = (0+$totalNbReads)/1000000;
    push(@nbReadsPerMillion, $readsPerMillion);
#    print $readsPerMillion."\n";
    
  }
  @samplesTotalNbReads = undef;
  
  while (<READS_FILE>) {
    chomp;
    my $line = $_;
    my @transcriptInfo = split(/\t/, $line);
    $transcriptInfo[0] =~ /.*_Length_(\d+)$/;
    my $transcriptLength = 0+$1;

    print $transcriptInfo[0];
    for(my $idx=1; $idx < @transcriptInfo; $idx++) {
      my $nbReads = 0+$transcriptInfo[$idx];
      printf "\t%.4f", ($nbReads/($nbReadsPerMillion[$idx-1] * $transcriptLength));
    }
    print "\n";
  }
  close(READS_FILE);
}
