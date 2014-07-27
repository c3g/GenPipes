#!/usr/bin/env perl

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

package TophatBowtie;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin";

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

  # Bowtie index basename
  my $bowtie_idx_basename = LoadConfig::getParam($rH_cfg, 'align', 'bowtieIndexBasename', 1, 'prefixpath'); ;
  
  # -G and --transcriptome-index options
  my $refFile      = LoadConfig::getParam($rH_cfg, 'align','referenceGtf',       0, 'filepath');
  my $refFileIndex = LoadConfig::getParam($rH_cfg, 'align','transcriptomeIndex', 0, 'prefixpath');
  my $refOption = ' ';
  if ($refFile) {
    $refOption .= '-G ' . $refFile;
	if ($refFileIndex) {
		$refOption .= ' --transcriptome-index ' . $refFileIndex;
	}
  }
  
  # Tophat any other options
  my $otherOptions      = LoadConfig::getParam($rH_cfg, 'align','tophatMoreOptions', 0);


  my $ro_job = new Job();
  if (defined($pair2)) {
    $ro_job->testInputOutputs([$pair1, $pair2], [$outputBAM]);
  } else {
    $ro_job->testInputOutputs([$pair1], [$outputBAM]);
  }

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [
      ['align', 'moduleVersion.bowtie'],
      ['align', 'moduleVersion.tophat'],
      ['align', 'moduleVersion.samtools']
    ]) . ' &&'; 
    $command .= ' tophat';
    $command .= ' --rg-library \"' . $rH_laneInfo->{'libraryBarcode'} . '\"';
    $command .= ' --rg-platform \"' . LoadConfig::getParam($rH_cfg, 'align', 'platform') . '\"';
    $command .= ' --rg-platform-unit \"' . $rH_laneInfo->{'lane'} . '\"';
    $command .= ' --rg-center \"' . LoadConfig::getParam($rH_cfg, 'align', 'TBInstitution') . '\"';
    $command .= ' --rg-sample ' . $sampleName;
    $command .= ' --rg-id ' . $rH_laneInfo->{'runId'};
    $command .= ' --library-type ' . LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
#     $command .= ' --fusion-search ' . LoadConfig::getParam($rH_cfg, 'align', 'fusionOption');
    $command .= ' -o ' . $laneDirectory;
    $command .= ' -p ' . LoadConfig::getParam($rH_cfg, 'align', 'TBAlnThreads', 1, 'int') . $refOption . ' ' . $otherOptions  ;
#     $command .= ' -g ' . LoadConfig::getParam($rH_cfg, 'align', 'maxReadLocation');

    $command .= ' ' . $bowtie_idx_basename;

    $command .= ' ' . $pair1;
    if (defined($pair2)) {
      $command .= ' ' . $pair2;
    }

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

1;
