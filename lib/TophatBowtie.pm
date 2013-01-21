#!/usr/env/perl

=head1 NAME

I<TophatBowtie>

=head1 SYNOPSIS

TophatBowtie-> align()

=head1 DESCRIPTION

B<TopHatBowtie> is a library to manage the different tools offerts by both tophat and Bowtie software

Input = file_name

Output = array


=head1 AUTHOR
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package TopHatBowtie;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------
sub align {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rH_laneInfo = shift;
  my $pair1       = shift;
  my $pair2       = shift;
  my $single1     = shift;
  my $single2     = shift;

  my $command = "";

  if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
    $command = singleCommand($rH_cfg, $sampleName, $rH_laneInfo, $single1);
  }
  elsif($rH_laneInfo->{'runType'} eq "PAIRED_END") {
    $command = pairCommand($rH_cfg, $sampleName, $rH_laneInfo, $pair1, $pair2);
  }
  else {
    die "Unknown runType: ".$rH_laneInfo->{' runType '}."\n";
  }
    
  return $command;
}

sub pairCommand {
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rH_laneInfo = shift;
  my $pair1       = shift;
  my $pair2       = shift;

  my $laneDirectory = "alignment/" . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
  my $outputBAM = $laneDirectory . 'accepted_hits.bam';
  my $bamFileDate = -M $outputBAM;

  my $commands;
  if (!defined($bamFileDate) || !defined(-M $pair1) || !defined(-M $pair2) || $bamFileDate < -M $pair1 || $bamFileDate < -M $pair2) {
    my $BTCommand = "";

    $BTCommand .= 'module load ' .LoadConfig::getParam($rH_cfg, 'align','bowtieModule') .' ;'; 
    $BTCommand .= '  module load ' .LoadConfig::getParam($rH_cfg, 'align','bowtieModule') .' ;'; 
    $BTCommand .= ' tophat';
    $BTCommand .= ' --rg-library \"' . $rH_laneInfo->{'libraryBarcode'} .'\"';
    $BTCommand .= ' --rg-platform \"' .LoadConfig::getParam($rH_cfg, 'align','platform') .'\"';
    $BTCommand .= ' --rg-platform-unit \"' .$rH_laneInfo->{'lane'} .'\"';
    $BTCommand .= ' --rg-center \"'. LoadConfig::getParam($rH_cfg, 'align','TBInstitution') .'\"';
    $BTCommand .= ' --rg-sample '. $sampleName;
    $BTCommand .= ' --rg-platform ' .$rH_laneInfo->{'runId'};
    $BTCommand .= ' --library-type '. LoadConfig::getParam($rH_cfg, 'align','strandInfo');
    $BTCommand .= ' --fusion-search '. LoadConfig::getParam($rH_cfg, 'align','fusionOption');
    $BTCommand .= ' -o ' .$laneDirectory;
    $BTCommand .= ' -p '. LoadConfig::getParam($rH_cfg, 'align','TBAlnThreads') .' -G '. LoadConfig::getParam($rH_cfg, 'align','referenceGtf');
    $BTCommand .= ' '. LoadConfig::getParam($rH_cfg, 'align','referenceFasta');
    $BTCommand .= ' '. $pair1 .' '. $pair2;
    $commands= $BTCommand;
  }

  return $commands;
}

sub singleCommand {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rH_laneInfo = shift;
  my $single      = shift;

  my $laneDirectory = "alignment/" . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
  my $outputBAM = $laneDirectory . 'accepted_hits.bam';
  my $bamFileDate = -M $outputBAM;

  my $commands;
  if (!defined($bamFileDate) || !defined(-M $pair1) || !defined(-M $pair2) || $bamFileDate < -M $pair1 || $bamFileDate < -M $pair2) {
    my $BTCommand = "";

    $BTCommand .= 'module load ' .LoadConfig::getParam($rH_cfg, 'align','bowtieModule') .' ;'; 
    $BTCommand .= '  module load ' .LoadConfig::getParam($rH_cfg, 'align','bowtieModule') .' ;'; 
    $BTCommand .= ' tophat';
    $BTCommand .= ' --rg-library \"' . $rH_laneInfo->{'libraryBarcode'} .'\"';
    $BTCommand .= ' --rg-platform \"' .LoadConfig::getParam($rH_cfg, 'align','platform') .'\"';
    $BTCommand .= ' --rg-platform-unit \"' .$rH_laneInfo->{'lane'} .'\"';
    $BTCommand .= ' --rg-center \"'. LoadConfig::getParam($rH_cfg, 'align','TBInstitution') .'\"';
    $BTCommand .= ' --rg-sample '. $sampleName;
    $BTCommand .= ' --rg-platform ' .$rH_laneInfo->{'runId'};
    $BTCommand .= ' --library-type '. LoadConfig::getParam($rH_cfg, 'align','strandInfo');
    $BTCommand .= ' --fusion-search '. LoadConfig::getParam($rH_cfg, 'align','fusionOption');
    $BTCommand .= ' -o ' .$laneDirectory;
    $BTCommand .= ' -p '. LoadConfig::getParam($rH_cfg, 'align','TBAlnThreads') .' -G '. LoadConfig::getParam($rH_cfg, 'align','referenceGtf');
    $BTCommand .= ' '. LoadConfig::getParam($rH_cfg, 'align','referenceFasta');
    $BTCommand .= ' '. $single;
    $commands= $BTCommand;
  }

  return $commands;
}
 

1;
