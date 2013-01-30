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

  my $laneDirectory = "alignment/" . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
  my $outputBAM = $laneDirectory . 'accepted_hits.bam';
  my $bamFileDate = -M $outputBAM;

  my $command;
  if (!defined($bamFileDate) || !defined(-M $ingle)  || $bamFileDate < -M $single) {
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'align','bowtieModule') .' ;'; 
    $command .= '  module load ' .LoadConfig::getParam($rH_cfg, 'align','bowtieModule') .' ;'; 
    $command .= ' tophat';
    $command .= ' --rg-library \"' . $rH_laneInfo->{'libraryBarcode'} .'\"';
    $command .= ' --rg-platform \"' .LoadConfig::getParam($rH_cfg, 'align','platform') .'\"';
    $command .= ' --rg-platform-unit \"' .$rH_laneInfo->{'lane'} .'\"';
    $command .= ' --rg-center \"'. LoadConfig::getParam($rH_cfg, 'align','TBInstitution') .'\"';
    $command .= ' --rg-sample '. $sampleName;
    $command .= ' --rg-platform ' .$rH_laneInfo->{'runId'};
    $command .= ' --library-type '. LoadConfig::getParam($rH_cfg, 'align','strandInfo');
    $command .= ' --fusion-search '. LoadConfig::getParam($rH_cfg, 'align','fusionOption');
    $command .= ' -o ' .$laneDirectory;
    $command .= ' -p '. LoadConfig::getParam($rH_cfg, 'align','TBAlnThreads') .' -G '. LoadConfig::getParam($rH_cfg, 'align','referenceGtf');
    $command .= ' '. LoadConfig::getParam($rH_cfg, 'align','referenceFasta');
    $command .= ' '. $pair1 .' ' .$pair2;
  }

  return $command;
}
 

1;
