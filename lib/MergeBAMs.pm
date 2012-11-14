#!/usr/env/perl

=head1 NAME

I<Trimmomatic>

=head1 SYNOPSIS

Trimmomatic->trim()

=head1 DESCRIPTION

B<Trimmomatic> is a library that trims fastqs

Input = file_name

Output = array


=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package MergeBAMs;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub merge {
  my $rH_cfg            = shift;
  my $sampleName        = shift;
  my $rAoH_sampleLanes  = shift;

  my $latestBam;
  my $bamInputs;
  my $countInputs;
  my $outputBAM = $sampleName.'/'.$sampleName.'.sorted.bam';
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $directory = $sampleName."/run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'}."/";
    my $sortedLaneBamFile = $directory.$rH_laneInfo->{'name'}.".sorted.bam";
    my $laneStatsFile = $directory.$rH_laneInfo->{'name'}.".counts";
    my $runName = $sampleName."_run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'};


    if(!defined($latestBam)) {
      $latestBam = -M $sortedLaneBamFile;
    }
    else {
      my $modDate = -M $sortedLaneBamFile;
      if($modDate > $latestBam) {
        $latestBam = $modDate;
      }
    }
    $bamInputs .= 'INPUT='.$sortedLaneBamFile.' ';
    $countInputs .= $laneStatsFile.' ';
  }

  my $command;
  if(!defined($latestBam) || !defined(-M $outputBAM) || $latestBam > -M $outputBAM) {
    $command .= 'module load mugqic/picard/1.77 ;';
    $command .= ' java -Xmx'.LoadConfig::getParam($rH_cfg, 'mergeLanes', 'mergeRam').' -jar ${PICARD_HOME}/MergeSamFiles.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true';
    $command .= ' '.$bamInputs;
    $command .= ' OUTPUT='.$outputBAM;
    $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'mergeLanes', 'mergeRecInRam');
    $command .= ' ; cat '.$countInputs.' > $sampleName/$sampleName.runLane.counts';
  }
  return $command;
}

1;
