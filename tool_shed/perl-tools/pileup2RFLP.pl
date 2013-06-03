#!/usr/bin/perl

#
# $Id$
#

use strict;
use warnings;
use Getopt::Std;
use Thread;
use Thread::Queue;
use Bio::Restriction::Analysis;
use Bio::Restriction::EnzymeCollection;
use Bio::DB::Fasta;
use Bio::Seq::SeqFactory;

use vars qw($VERSION);
$VERSION = '1.0';

my $pileupQueue = new Thread::Queue;
my $enzymeQueue = new Thread::Queue;
my $factory     = Bio::Seq::SeqFactory->new();
my $flankSize = 1000;
my $rebase_collection;
my $output;

my %companies = (
  'B' => 'Invitrogen Corporation (4/10)',
  'C' => 'Minotech Biotechnology (4/10)',
  'E' => 'Stratagene (3/10)',
  'F' => 'Fermentas International Inc. (4/10)',
  'H' => 'American Allied Biochemical, Inc. (4/10)',
  'I' => 'SibEnzyme Ltd. (4/10)',
  'J' => 'Nippon Gene Co., Ltd. (4/10)',
  'K' => 'Takara Bio Inc. (4/10)',
  'M' => 'Roche Applied Science (4/10)',
  'N' => 'New England Biolabs (4/10)',
  'O' => 'Toyobo Biochemicals (9/09)',
  'Q' => 'Molecular Biology Resources - CHIMERx (4/10)',
  'R' => 'Promega Corporation (3/10)',
  'S' => 'Sigma Chemical Corporation (4/10)',
  'U' => 'Bangalore Genei (4/10)',
  'V' => 'Vivantis Technologies (4/10)',
  'X' => 'EURx Ltd. (4/10)',
  'Y' => 'CinnaGen Inc. (4/10)'
);

&main();

sub main() {
  my %opts = ( r => undef, e => undef, p => undef, t => 1, o => 'output.out' );
  getopts( 'r:e:p:o:t:', \%opts );

  die(
    qq/
Usage:   pileup2rflp.pl [options] -r <fasta_reference> [-e <rebase_enzymes>] [-t nbThreads] -p <pileup>
Version: $VERSION
Options: -r File   Reference sequences. Name of the sequences must match those found in the pileup file.
         -p File   SNPs in pileup format.
         -t INT    Number of threads to use [$opts{t}].
         -e File   ReBase enzyme file.
         -o File   Output file [$opts{o}]
\n/
  ) if ( !defined( $opts{r} ) || !defined( $opts{p} ) );
  if ( $opts{t} < 1 ) {
    die "Nb Threads must be > 0. You set $opts{t}";
  }

  if($opts{o} eq '-') {
    $output = \*STDOUT;
  }
  else {
    open($output, ">$opts{o}") or die "Couldn't write to file '$opts{o}'\n";
  }

  if ( !defined( $opts{e} ) ) {
    $rebase_collection = Bio::Restriction::EnzymeCollection->new;
  }
  else {
    print "Loading enzymes from file $opts{e}\n";
    my $rebase = Bio::Restriction::IO->new(
      -file   => $opts{e},
      -format => 'withrefm'
    );
    $rebase_collection = $rebase->read();
  }

  open( PILEUP, "<$opts{p}" ) or die "Couldn't open pileup file '$opts{p}'\n";
  my $total = 0;
  while ( my $line = <PILEUP> ) {
    $total++;
    chomp($line);

    my ( $chr, $position, $refAllele, $snpCode, @values ) = split( /\t/, $line );
    my %pileupRow = (
      'chr'       => $chr,
      'position'  => $position,
      'refAllele' => $refAllele,
      'snpCode'   => $snpCode,
      'other'     => \@values
    );
    $pileupQueue->enqueue( \%pileupRow );
  }
  print "Total SNPs: $total\n";
  
  my @threads;
  for(my $i=0; $i < $opts{t}; $i++) {
    push(@threads, new Thread(\&analyseSNP, $opts{r}));
  }

  my $i=0;
  while(scalar(@threads) > 0) {
    print "Waiting for thread $i\n";
    my $thread = shift(@threads);
    $thread->join();
    print "Thread $i done\n";
    $i++;
  }
  
#  analyseSNP();
#  print "\n";
  
  outputWriter();
  
  if($opts{o} ne '-') {
    close($output);
  }
  print "Finished!\n";
}

sub getVendorName {
  my $enzyme = shift;
  my @vendorNames;

  my @vendorsList = $enzyme->vendors();
  if(scalar(@vendorsList) == 1) {
    my @vendors = @{$vendorsList[0]};
    foreach my $vendorKey (@vendors) {
      if(!defined($companies{$vendorKey})) {
        die "Unknown vendor: $vendorKey\n";
      }
      push(@vendorNames, $companies{$vendorKey});
    }
  }
  elsif (scalar(@vendorsList) > 1) {
    die "More than 1 vendor found: ".$enzyme->name()."\n";
  }
  return join(",", @vendorNames);
}

sub outputWriter {
  while(my $rH_enzyme = $enzymeQueue->dequeue_nb()) {
    if(!defined($rH_enzyme)) {
      last;
    }

    printf $output $rH_enzyme->{'chr'};
    printf $output "\t";
    printf $output $rH_enzyme->{'position'};
    printf $output "\t";
    printf $output $rH_enzyme->{'enzymeName'};
    printf $output "\t";
    printf $output $rH_enzyme->{'ref'};
    printf $output "\t";
    printf $output $rH_enzyme->{'mutant'};
    printf $output "\t";
    printf $output $rH_enzyme->{'vendors'};
    printf $output "\n";
  }  
}

sub analyseSNP {
  my $refDbPath = shift;
  my $seq = undef;
  # create database from directory of fasta files
  my $refDb = Bio::DB::Fasta->new( $refDbPath );

  while(my $rH_pileupRow  = $pileupQueue->dequeue_nb()) {
    if(!defined($rH_pileupRow)) {
      last;
    }
    
    my $chr           = $rH_pileupRow->{'chr'};
    my $position      = $rH_pileupRow->{'position'};
    my $refAllele     = $rH_pileupRow->{'refAllele'};
    my $snpCode       = $rH_pileupRow->{'snpCode'};

    print "Remaining: ".$pileupQueue->pending()."\n";

    if ( $refAllele eq '*' ) {
      #warn "We don't support indels";
      next;
    }
  
    if ( !defined($seq) || $seq->display_id() ne $chr ) {
      $seq = $refDb->get_Seq_by_id($chr);
      if ( !defined($seq) ) {
        my @ids = $refDb->get_all_ids();
  
        #warn "Couldn't find $chr. Valid ids: @ids\n";
        next;
      }
      print "Changed Chromosome to ".$chr."\n";
    }

    print "Converting IUPAC to SNP allele\n";
    my %IUB     = Bio::Tools::IUPAC->iupac_iub();
    my @alleles = @{ $IUB{$snpCode} };
    if ( scalar(@alleles) > 2 ) {
      die "We don't support tri or quadri allelic snps\n";
    }
  
    my $snpAllele;
    if ( scalar(@alleles) == 1 ) {
      $snpAllele = $alleles[0];
    }
    else {
      if ( $alleles[0] eq uc($refAllele) ) {
        $snpAllele = $alleles[1];
      }
      else {
        $snpAllele = $alleles[0];
      }
    }
  
    print "Create sequence flank5+SNP+flank3\n";
  
    my $startFlank = ( $position - $flankSize ) < 0 ? 0 : $position - $flankSize;
    my $stopFlank = ( $position + $flankSize ) > $seq->length() ? $seq->length() : $position + $flankSize;
  
    my $flankedSeq = $seq->subseq( $startFlank, $stopFlank );
    my $flank5 = substr( $flankedSeq, 0, $flankSize );
    my $flank3 = substr( $flankedSeq, $flankSize + 1 );    # +1 to remove the snp
    my $flankedMutantSeq = $flank5 . $snpAllele . $flank3;
  
    # Debug flanks
    #print "Values:\n";
    #print "$flankedSeq\n";
    #print "$flankedMutantSeq\n";
    #print $refAllele."\n";
    #print "@alleles\n";
    #print $flank5." ".$snpAllele. " ".$flank3."\n";
  
    my $origSeq = $factory->create(
      -seq => $flankedSeq,
      -id  => 'orig'
    );
    my $mutaSeq = $factory->create(
      -seq => $flankedMutantSeq,
      -id  => 'mutant'
    );
  
    print "Analyse reference\n";
    my $origAnalysis = Bio::Restriction::Analysis->new(
      -seq     => $origSeq,
      -enzymes => $rebase_collection
    );
    print "Analyse mutant\n";
    my $mutaAnalysis = Bio::Restriction::Analysis->new(
      -seq     => $mutaSeq,
      -enzymes => $rebase_collection
    );
  
    print "Testing enzymes";
    my $allCutters = $origAnalysis->cutters();
    
    my $allMutantCutters = $mutaAnalysis->cutters();
    my %mutantCutterEnzymes;
    foreach my $enz ( $allMutantCutters->each_enzyme ) {
      $mutantCutterEnzymes{$enz->name()} = 1;
    }
    
    foreach my $enz ( $allCutters->each_enzyme ) {
      print ".";
      
      delete $mutantCutterEnzymes{$enz->name()};
      my @vendors = $enz->vendors();
      my $vendorNames = getVendorName($enz);
      
      if(length($vendorNames) == 0) {
        next;
      }

      my $origFragSize;
      my $mutaFragSize;
      eval {
        $origFragSize = $origAnalysis->fragments($enz);
        $mutaFragSize = $mutaAnalysis->fragments($enz);
      };
      if( $@ ) {
        print STDERR "Caught exception\n";
        next;
      }
  
      if ( $origFragSize != $mutaFragSize ) {
        my %enzyme = (
          'chr' => $chr,
          'position' => $position,
          'enzymeName' => $enz->name(),
          'ref' => $origFragSize,
          'mutant' => $mutaFragSize,
          'vendors' => getVendorName($enz)
        );
        $enzymeQueue->enqueue(\%enzyme);
      }
    }
    
    # If the mutant has some cutters the original doesn't, use them
    foreach my $enzymeName ( keys(%mutantCutterEnzymes) ) {
      print ",";
      my $enz = $allMutantCutters->get_enzyme($enzymeName);
      
      my @vendors = $enz->vendors();
      my $vendorNames = getVendorName($enz);
      
      if(length($vendorNames) == 0) {
        next;
      }

      my $origFragSize;
      my $mutaFragSize;
      eval {
        $origFragSize = $origAnalysis->fragments($enz);
        $mutaFragSize = $mutaAnalysis->fragments($enz);
      };
      if( $@ ) {
        print STDERR "Caught exception\n";
        next;
      }

      my %enzyme = (
        'chr' => $chr,
        'position' => $position,
        'enzymeName' => $enz->name(),
        'ref' => $origFragSize,
        'mutant' => $mutaFragSize,
        'vendors' => getVendorName($enz)
      );
      $enzymeQueue->enqueue(\%enzyme);
    }    
    print "Done\n";
    
  }
}

