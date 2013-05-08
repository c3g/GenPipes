#!/usr/bin/perl -w

use strict;
use Set::IntervalTree;

&main();

sub main {
  my %peaks;
#  my $repeats = Set::IntervalTree->new;
#  my $genes = Set::IntervalTree->new;
  my %repeats;
  my %genes;
  parsePeaks( \%peaks );
  getRepeats( \%repeats );
  getGenes( \%genes );

  my $removedForTags=0;
  my $removedForRpt=0;
  my $removedForGeneDistance=0;
  my $macsPeakId=1;
  for my $rH_peaks ( values %peaks ) {
    print STDERR "Testing peak: ".$rH_peaks->{'chr'}.":".$rH_peaks->{'start'}."\n";
    my $repeatsFound = $repeats{$rH_peaks->{'chr'}}->fetch($rH_peaks->{'start'}, $rH_peaks->{'end'});
    my $genesFound = $genes{$rH_peaks->{'chr'}}->fetch($rH_peaks->{'start'}, $rH_peaks->{'end'});
    
    my $removed=0;
    if($rH_peaks->{'tags'} < 4) {
      $removedForTags++;
      print STDERR "Tags removed: ".$rH_peaks->{'chr'}.":".$rH_peaks->{'start'}."\n";
      $removed=1;
    }
    
    if(@$repeatsFound > 0) {
      my $midPeakLength = ($rH_peaks->{'end'} - $rH_peaks->{'start'})/2;
      for my $rH_repeatFound (@$repeatsFound) {
        if($removed==1) {
          last;
        }
        my $midRptLength = ($rH_repeatFound->{'end'} - $rH_repeatFound->{'start'})/2;
        my $overlapStart = $rH_peaks->{'start'};
        if($rH_peaks->{'start'} < $rH_repeatFound->{'start'}) {
          $overlapStart = $rH_repeatFound->{'start'};
        }
        my $overlapEnd = $rH_peaks->{'end'};
        if($rH_peaks->{'end'} > $rH_repeatFound->{'end'}) {
          $overlapEnd = $rH_repeatFound->{'end'};
        }

				my $overlapLength = $overlapEnd - $overlapStart;

        if($overlapLength > $midRptLength || $overlapLength > $midPeakLength) {
          $removedForRpt++;
          print STDERR "Rpt removed: ".$rH_peaks->{'chr'}.":".$rH_peaks->{'start'}."\n";
          $removed=1;
        }
#        if($rH_peaks->{'start'} > $rH_repeatFound->{'start'}) {
#          my $overlapLength = $rH_repeatFound->{'end'} - $rH_peaks->{'start'};
#          if($overlapLength >= $rptLength) {
#            $removedForRpt++;
#            $removed=1;
#          }
#        }
#        else {
#          my $overlapLength = $rH_peaks->{'end'} - $rH_repeatFound->{'start'};
#          if($overlapLength >= $rptLength) {
#            $removedForRpt++;
#            $removed=1;
#          }
#        }
      }
    }
    
    if($removed == 0 && @$genesFound > 0) {
      if($$genesFound[0]->{'strand'} eq '+') {
        if($rH_peaks->{'end'} < $$genesFound[0]->{'start'}-10000 || $rH_peaks->{'start'} > $$genesFound[0]->{'start'}+1000) {
          $removedForGeneDistance++;
          print STDERR "Gene removed: ".$rH_peaks->{'chr'}.":".$rH_peaks->{'start'}."\n";
          $removed=1;
        }
      }
      else {
        if($rH_peaks->{'start'} > $$genesFound[0]->{'end'}+10000 || $rH_peaks->{'end'} < $$genesFound[0]->{'end'}-1000) {
          $removedForGeneDistance++;
          print STDERR "Gene removed: ".$rH_peaks->{'chr'}.":".$rH_peaks->{'start'}."\n";
          $removed=1;
        }
      }
    }
    
    if($removed == 0) {
      print "xls\t";
      print $rH_peaks->{'chr'};
      print "\t";
      print $rH_peaks->{'start'};
      print "\t";
      print $rH_peaks->{'end'};
      print "\t";
      print $rH_peaks->{'length'};
      print "\t";
      print $rH_peaks->{'summit'};
      print "\t";
      print $rH_peaks->{'tags'};
      print "\t";
      print $rH_peaks->{'pValue'};
      print "\t";
      print $rH_peaks->{'foldEnrichment'};
      print "\n";
      print "bed\t";
      print $rH_peaks->{'chr'};
      print "\t";
      if(($rH_peaks->{'start'}-1) < 0) {
        print "0";
      }
      else {
        print $rH_peaks->{'start'}-1;
      }
      print "\t";
      print $rH_peaks->{'end'};
      print "\tMACS_peak_";
      print $macsPeakId;
      print "\t";
      print $rH_peaks->{'pValue'};
      print "\n";
			
    }
    $macsPeakId++;
  }
  print STDERR "Nb removed for nb tags:                           ".$removedForTags."\n";
  print STDERR "Nb removed for covering 50% of a repeat:          ".$removedForRpt."\n";
  print STDERR "Nb removed for being outside gene start - 2500bp: ".$removedForGeneDistance."\n";
}

sub parsePeaks() {
  my $rHoH_peaks = shift;

  my $macsOutput = $ARGV[0];
  open( MACS, "<$macsOutput" ) or die "Can't open macs file\n";
  while (<MACS>) {
    chomp;
    if ( substr( $_, 0, 1 ) eq '#' || length($_) == 0) {
      next;
    }

		if(/^chr\tstart/) {
			next;
		}

    my @results = split;
    my %peak;
    $peak{'chr'}            = $results[0];
    $peak{'start'}          = 0+$results[1];
    $peak{'end'}            = 0+$results[2];
    $peak{'length'}         = $results[3];
    $peak{'summit'}         = $results[4];
    $peak{'tags'}           = $results[5];
    $peak{'pValue'}         = $results[6];
    $peak{'foldEnrichment'} = $results[7];
    print STDERR "Read peak: ".$peak{'chr'}.":".$peak{'start'}."\n";
    $rHoH_peaks->{ $peak{'chr'}.'-'.$peak{'start'} } = \%peak;
  }
  print STDERR "Read ".keys(%$rHoH_peaks)." peaks\n";
  close(MACS);
}

sub getRepeats {
  my $rHoO_repeats = shift;

  my $repeatMsk = $ARGV[1];
  open( RMSK, "<$repeatMsk" ) or die "Can't open repeatMsk file\n";
  while (<RMSK>) {
    chomp;

    my @results = split;
    my %rmsk;
    $rmsk{'chr'}            = $results[5];
    $rmsk{'start'}          = 0+$results[6];
    $rmsk{'end'}            = 0+$results[7];
    $rmsk{'strand'}         = $results[9];
    $rmsk{'name'}         = $results[10];

    if(!defined $rHoO_repeats->{ $rmsk{'chr'} }) {
      $rHoO_repeats->{ $rmsk{'chr'} } = Set::IntervalTree->new;
    }
    $rHoO_repeats->{ $rmsk{'chr'} }->insert(\%rmsk,$rmsk{'start'}+1,$rmsk{'end'});
  }
  close(RMSK);

  $repeatMsk = $ARGV[2];
  open( RMSK, "<$repeatMsk" ) or die "Can't open TRF file\n";
  while (<RMSK>) {
    chomp;

    my @results = split;
    my %rmsk;
    $rmsk{'chr'}            = $results[1];
    $rmsk{'start'}          = 0+$results[2];
    $rmsk{'end'}            = 0+$results[3];
    $rmsk{'name'}         = $results[4];
    if(!defined $rHoO_repeats->{ $rmsk{'chr'} }) {
      $rHoO_repeats->{ $rmsk{'chr'} } = Set::IntervalTree->new;
    }
    $rHoO_repeats->{ $rmsk{'chr'} }->insert(\%rmsk,$rmsk{'start'}+1,$rmsk{'end'});
  }
  close(RMSK);
}

sub getGenes {
  my $rHoO_genes = shift;

  my $gene = $ARGV[3];
  open( GENE, "<$gene" ) or die "Can't open refGene file\n";
  while (<GENE>) {
    chomp;

    my @results = split;
    my %gene;
    $gene{'chr'}            = $results[2];
    $gene{'strand'}         = $results[3];
    $gene{'start'}          = 0+$results[4];
    $gene{'end'}            = 0+$results[5];
    $gene{'name'}         = $results[12];
    
    if(!defined $rHoO_genes->{ $gene{'chr'} }) {
      $rHoO_genes->{ $gene{'chr'} } = Set::IntervalTree->new;
    }
    $rHoO_genes->{ $gene{'chr'} }->insert(\%gene,$gene{'start'}+1,$gene{'end'});
    
  }
  close(GENE);
}

__END__
