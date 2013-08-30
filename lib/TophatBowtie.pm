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

package TophatBowtie;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use LoadConfig;
use PipelineUtils;

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

#   #------ flefebvr Tue 16 Apr 09:04:33 2013 
#   # Bowtie index basename (assumes reference fasta is named basename.extension) 
#   (my $bwa_idx_basename  = LoadConfig::getParam($rH_cfg, 'align','bowtieRefIndex') ) =~ s/\.[^.]+$//;
#   #------
  ####mbourgey - Francois' change does not works when using hg1k

  my $bwa_idx_basename ;
  if (-e LoadConfig::getParam($rH_cfg, 'align','referenceFasta') .'.1.bt2') {
     $bwa_idx_basename  = LoadConfig::getParam($rH_cfg, 'align','referenceFasta');
  } else {
     ($bwa_idx_basename  = LoadConfig::getParam($rH_cfg, 'align','referenceFasta') ) =~ s/\.[^.]+$//;
  }
  my $refFile=LoadConfig::getParam($rH_cfg, 'align','referenceGtf');
  my $refOption=' ';
  if ($refFile ne ' ') {
    $refOption .= '-G '.$refFile;
  }

  my $up2date = PipelineUtils::testInputOutputs([$pair1, $pair2], [$outputBAM]);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'align','moduleVersion.bowtie') ; 
    $command .= ' ' .LoadConfig::getParam($rH_cfg, 'align','moduleVersion.tophat') ;
    $command .= ' ' .LoadConfig::getParam($rH_cfg, 'align','moduleVersion.samtools').' &&'; 
    $command .= ' tophat';
    $command .= ' --rg-library \"' . $rH_laneInfo->{'libraryBarcode'} .'\"';
    $command .= ' --rg-platform \"' .LoadConfig::getParam($rH_cfg, 'align','platform') .'\"';
    $command .= ' --rg-platform-unit \"' .$rH_laneInfo->{'lane'} .'\"';
    $command .= ' --rg-center \"'. LoadConfig::getParam($rH_cfg, 'align','TBInstitution') .'\"';
    $command .= ' --rg-sample '. $sampleName;
    $command .= ' --rg-id ' .$rH_laneInfo->{'runId'};
    $command .= ' --library-type '. LoadConfig::getParam($rH_cfg, 'align','strandInfo');
#     $command .= ' --fusion-search '. LoadConfig::getParam($rH_cfg, 'align','fusionOption');
    $command .= ' -o ' .$laneDirectory;
    $command .= ' -p '. LoadConfig::getParam($rH_cfg, 'align','TBAlnThreads') .$refOption;
#     $command .= ' -g '. LoadConfig::getParam($rH_cfg, 'align','maxReadLocation');

    #------ flefebvr Tue 16 Apr 09:04:54 2013 
    #$command .= ' '. LoadConfig::getParam($rH_cfg, 'align','bowtieRefIndex');
    $command .= ' '. $bwa_idx_basename;
    #------

    $command .= ' '. $pair1 .' ' .$pair2;
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

1;
