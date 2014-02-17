#!/usr/env/perl

=head1 NAME

I<SVTools>

=head1 SYNOPSIS

Picard->merge()

=head1 DESCRIPTION

B<SVTools> is a library to analyse BAMs for Structural Variants

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package SVtools;

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
sub runPairedDNAC {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBins     = shift;
  my $outputPrefix  = shift;
  my $window        = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputBins], [$outputPrefix . '.txt']);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['runPairedDNAC', 'moduleVersion.cranR'], ['runPairedDNAC', 'moduleVersion.svtools']]) . ' &&';
    $command .= ' Rscript \${SVTOOLS_HOME}/Cancer/RunDNAC.6.0.R';
    $command .= ' -f ' . $inputBins;
    $command .= ' -b ' . $window;
    $command .= ' -o ' . $outputPrefix;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub filterDNAC {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $inputDNACCalls  = shift;
  my $outputPrefix    = shift;
  my $cnvProx         = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputDNACCalls], [$outputPrefix . '.filteredSV.txt']);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['filterSV', 'moduleVersion.svtools']]) . ' &&';
    $command .= ' \${SVTOOLS_HOME}/Cancer/filterOutDNAC.sh';
    $command .= ' ' . $inputDNACCalls;
    $command .= ' ' . $outputPrefix . '.txt';
    $command .= ' ' . $sampleName;
    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'filterSV', 'minBinCNV', 1, 'int');
    $command .= ' && ';
    $command .= filterResults($rH_cfg, 'filterSV', $outputPrefix, $cnvProx);
    $command .= ' && ';
    $command .= generateBedResults($rH_cfg, $outputPrefix);

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub filterBrD {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBrDCalls = shift;
  my $outputPrefix  = shift;
  my $normalFile    = shift;
  my $tumorFile     = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputBrDCalls], [$outputPrefix . '.filteredSV.txt']);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['filterSVC', 'moduleVersion.svtools']]) . ' &&';
    $command .= ' \${SVTOOLS_HOME}/Cancer/filterOutBrD.py';
    $command .= ' -f ' . $inputBrDCalls;
    $command .= ' -o ' . $outputPrefix . '.txt';
    $command .= ' -s ' . $sampleName;
    $command .= ' -n ' . LoadConfig::getParam($rH_cfg, 'filterSV', 'minReadSupport', 1, 'int');
    $command .= ' -t ' . LoadConfig::getParam($rH_cfg, 'filterSV', 'minReadSupport', 1, 'int');
    $command .= ' -b ' . $normalFile;
    $command .= ' -c ' . $tumorFile;
    $command .= ' && ';
    $command .= filterResults($rH_cfg, 'filterSV', $outputPrefix, '');
    $command .= ' && ';
    $command .= generateBedResults($rH_cfg, $outputPrefix);

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub filterPI {
  my $rH_cfg       = shift;
  my $sampleName   = shift;
  my $outputPrefix = shift;
  my $normalFile   = shift;
  my $tumorFile    = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$outputPrefix], [$outputPrefix . '.filteredSV.txt']);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['filterSV', 'moduleVersion.svtools']]) . ' &&';
    $command .= ' \${SVTOOLS_HOME}/Cancer/filterOutPI.py';
    $command .= ' -f ' . $outputPrefix;
    $command .= ' -o ' . $outputPrefix . '.txt';
    $command .= ' -s ' . $sampleName;
    $command .= ' -n ' . LoadConfig::getParam($rH_cfg, 'filterSV', 'minReadSupport', 1, 'int');
    $command .= ' -t ' . LoadConfig::getParam($rH_cfg, 'filterSV', 'minReadSupport', 1, 'int');
    $command .= ' && ';
    $command .= filterResults($rH_cfg, 'filterSV', $outputPrefix, '');
    $command .= ' && ';
    $command .= generateBedResults($rH_cfg, $outputPrefix);

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub filterResults {
  my $rH_cfg        = shift;
  my $stepIniPrefix = shift;
  my $outputPrefix  = shift;
  my $cnvProx       = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$outputPrefix . '.bed'], [$outputPrefix . '.txt']);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [[$stepIniPrefix, 'moduleVersion.svtools']]) . ' &&';
    $command .= ' \${SVTOOLS_HOME}/Cancer/filterBedResults.sh';
    $command .= ' ' . $outputPrefix . '.txt';
    $command .= ' ' . LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceMappabilityBed', 1, 'filepath');
    $command .= ' ' . LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceGeneCoordinates', 1, 'filepath');
    $command .= ' ' . LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceDGVCoordinates', 1, 'filepath');
    $command .= ' ' . LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceMicrosatellitesCoordinates', 1, 'filepath');
    $command .= ' ' . LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceRepeatMaskerCoordinates', 1, 'filepath');
    $command .= ' ' . LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceGenomeLengths', 1, 'filepath');
    $command .= ' ' . $outputPrefix . '.bed';
    $command .= ' ' . $outputPrefix . '.tmp';
    $command .= ' ' . $cnvProx;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub generateBedResults {
  my $rH_cfg          = shift;
  my $outputPrefix    = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$outputPrefix . '.bed.other.filteredSV.annotate.txt', $outputPrefix . '.bed.TumS.filteredSV.annotate.txt'], [$outputPrefix . '.bed.other.filteredSV.annotate.bed', $outputPrefix . '.bed.TumS.filteredSV.annotate.bed']);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['bedSV', 'moduleVersion.svtools']]) . ' &&';
    $command .= ' \${SVTOOLS_HOME}/Cancer/rtxt2rbed.sh ';
    $command .= ' ' . $outputPrefix . '.bed.other.filteredSV.annotate.txt';
    $command .= ' ' . $outputPrefix . '.bed.other.filteredSV.annotate.bed';
    $command .= ' && \${SVTOOLS_HOME}/Cancer/rtxt2rbed.sh';
    $command .= ' ' . $outputPrefix . '.bed.TumS.filteredSV.annotate.txt';
    $command .= ' ' . $outputPrefix . '.bed.TumS.filteredSV.annotate.bed';

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

1;
