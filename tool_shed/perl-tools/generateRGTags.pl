#!/usr/bin/perl

use strict;

&main();

sub main {
  my $samplesDescriptionFile = $ARGV[0];

  my $rA_SampleInfos = parseSheet($samplesDescriptionFile);

  
  for my $rH_Sample (@$rA_SampleInfos) {

    my $directory = $rH_Sample->{'name'}."/run".$rH_Sample->{'runId'}."_".$rH_Sample->{'lane'}."/";
    my $rgId = $rH_Sample->{'libraryBarcode'}."_".$rH_Sample->{'runId'}."_".$rH_Sample->{'lane'};
    my $rgTag = '@RG\tID:'.$rgId.'\tSM:'.$rH_Sample->{'name'}.'\tLB:'.$rH_Sample->{'libraryBarcode'}.'\tPU:run'.$rH_Sample->{'runId'}."_".$rH_Sample->{'lane'}.'\tCN:McGill University and Genome Quebec Innovation Center\tPL:Illumina';

    print $rH_Sample->{'name'}."\t".$rgTag."\n";
  }
}

sub parseSheet {
  my $fileName = shift;

  my @retVal;
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";
  <SAMPLE_SHEET>; # Ignore header
  while(<SAMPLE_SHEET>) {
    chomp;
    s/"//g;
    my @values = split(/,/, $_);

    my %sampleInfo;
    $sampleInfo{'name'} = $values[0];
    $sampleInfo{'libraryBarcode'} = $values[3];
    $sampleInfo{'runId'} = 0+$values[9];
    $sampleInfo{'lane'} = 0+$values[10];
    push(@retVal, \%sampleInfo);
  }

  return \@retVal;
}

