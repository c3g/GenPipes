#!/usr/bin/perl

use strict;

&main();

sub main {

  my %transcriptCount;
  my $nbSamplesUsed = -1;
  open(COUNT, $ARGV[0]) or die "Couldn't open ".$ARGV[0]."\n";
  while(<COUNT>){
    chomp;
    if(length($_) == 0) {
      next;
    }
    my @values = split(/\t/, $_);
    
    my @counts;
    my $idx=1;
    for(; $idx < @values; $idx++) {
      push(@counts, $values[$idx]);
    }
    $nbSamplesUsed = $idx-1;
    $transcriptCount{$values[0]} = \@counts;
  }
  close(COUNT);


  open(BLAST, $ARGV[1]) or die "Couldn't open ".$ARGV[1]."\n";
  while(<BLAST>){
    chomp;
    my @values = split(/\t/, $_);
    my $rA_count = $transcriptCount{$values[0]};

    print $values[0];    
    if(defined $rA_count) {
      for my $count (@{$rA_count}) {
        print "\t";    
        print $count;
      }
      delete $transcriptCount{$values[0]};
    }
    else {
      for(my $idx=0; $idx < $nbSamplesUsed; $idx++) {
        print "\tNA";
      }
    }

    for(my $idx=1; $idx < @values; $idx++) {
        print "\t";    
        print $values[$idx];
    }
    print "\n";
  }
  close(BLAST);

  for my $remainingTranscriptName (keys %transcriptCount) {
    print $remainingTranscriptName;
    for my $count (@{$transcriptCount{$remainingTranscriptName}}) {
        print "\t";
        print $count;
    }
    print "\n";
  }
}
