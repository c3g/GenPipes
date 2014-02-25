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

package Trimmomatic;

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
use Job;
use LoadConfig;

# SUB
#-----------------------

use constant {
  PAIR1_OUTPUT => 'pair1',
  PAIR2_OUTPUT => 'pair2',
  SINGLE1_OUTPUT => 'single1',
  SINGLE2_OUTPUT => 'single2',
};

sub trim {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rH_laneInfo = shift;
  my $outputDir   = shift;

  my $ro_job;

  my $rawReadDir = LoadConfig::getParam($rH_cfg, 'trim', 'rawReadDir', 1, 'dirpath');
  my $inputFastqPair1Name = $rawReadDir . '/' . $sampleName . '/run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '/' . $rH_laneInfo->{'read1File'};
  my $inputFastqPair2Name = undef;

  my $skipTrimming = LoadConfig::getParam($rH_cfg, 'trim', 'skip', 0);

  if ($rH_laneInfo->{'runType'} eq "SINGLE_END") {
    if (defined($skipTrimming) && $skipTrimming eq '1') {
      $ro_job = new Job();
      $ro_job->setOutputFileHash({SINGLE1_OUTPUT => $inputFastqPair1Name});
      $ro_job->setUp2Date(1);
    } else {
      $ro_job = singleCommand($rH_cfg, $sampleName, $rH_laneInfo, $outputDir);
    }
  } elsif ($rH_laneInfo->{'runType'} eq "PAIRED_END") {
    if (defined($skipTrimming) && $skipTrimming eq '1') {
      $inputFastqPair2Name = $rawReadDir . '/' . $sampleName . '/run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '/' . $rH_laneInfo->{'read2File'};
      $ro_job = new Job();
      $ro_job->setOutputFileHash({SINGLE1_OUTPUT => $inputFastqPair1Name, PAIR1_OUTPUT => $inputFastqPair1Name, PAIR2_OUTPUT => $inputFastqPair2Name});
      $ro_job->setUp2Date(1);
    } else {
      $ro_job = pairCommand($rH_cfg, $sampleName, $rH_laneInfo, $outputDir);
    }
  } else {
    die "Unknown runType: " . $rH_laneInfo->{'runType'} . "\n";
  }

  return $ro_job;
}

sub pairCommand {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rH_laneInfo = shift;
  my $outputDir   = shift;

  my $minQuality  = LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int');
  my $minLength   = LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int');
  my $adapterFile = LoadConfig::getParam($rH_cfg, 'trim', 'adapterFile', 1, 'filepath');
  my $headcrop = LoadConfig::getParam($rH_cfg, 'trim', 'headcrop', 0, 'int');

  my $rawReadDir = LoadConfig::getParam($rH_cfg, 'trim', 'rawReadDir', 1, 'dirpath');

  my $outputFastqPair1Name = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . $minQuality . 'l' . $minLength . '.pair1.fastq.gz';
  my $outputFastqPair2Name = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . $minQuality . 'l' . $minLength . '.pair2.fastq.gz';
  my $outputFastqSingle1Name = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . $minQuality . 'l' . $minLength . '.single1.fastq.gz';
  my $outputFastqSingle2Name = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . $minQuality . 'l' . $minLength . '.single2.fastq.gz';
  my $outputTrimLog = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.trim.out';
  my $outputTrimStats = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.trim.stats.csv';

  my $inputFastqPair1Name = $rawReadDir . '/' . $sampleName . '/run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '/' . $rH_laneInfo->{'read1File'};
  my $inputFastqPair2Name = $rawReadDir . '/' . $sampleName . '/run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '/' . $rH_laneInfo->{'read2File'};

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputFastqPair1Name, $inputFastqPair2Name], [$outputFastqPair1Name, $outputFastqPair2Name, $outputFastqSingle1Name, $outputFastqSingle2Name, $outputTrimLog, $outputTrimStats]);

  $ro_job->setOutputFileHash({PAIR1_OUTPUT => $outputFastqPair1Name, PAIR2_OUTPUT => $outputFastqPair2Name, SINGLE1_OUTPUT => $outputFastqSingle1Name, SINGLE2_OUTPUT => $outputFastqSingle2Name});

  if (!$ro_job->isUp2Date()) {
    my $command = "";
    $command .= LoadConfig::moduleLoad($rH_cfg, [
      ['trim', 'moduleVersion.java'],
      ['trim', 'moduleVersion.trimmomatic']
    ]);
    $command .= ' && java -XX:ParallelGCThreads=1 -Xmx2G -cp \$TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE';
    $command .= ' -threads ' . LoadConfig::getParam($rH_cfg, 'trim', 'nbThreads', 1, 'int');
    if ($rH_laneInfo->{'qualOffset'} eq "64") {
      $command .= ' -phred64';
    } else {
      $command .= ' -phred33';
    }
    $command .= ' ' . $inputFastqPair1Name;
    $command .= ' ' . $inputFastqPair2Name;
    $command .= ' ' . $outputFastqPair1Name . ' ' . $outputFastqSingle1Name;
    $command .= ' ' . $outputFastqPair2Name . ' ' . $outputFastqSingle2Name;
    if ($rH_laneInfo->{'qualOffset'} eq "64") {
      $command .= ' TOPHRED33';
    }
    if (defined($headcrop) && length($headcrop) > 0 && $headcrop > 0) {
      $command .= ' HEADCROP:' . $headcrop;
    }
    $command .= ' ILLUMINACLIP:' . $adapterFile . LoadConfig::getParam($rH_cfg, 'trim', 'clipSettings');
    if ($minQuality > 0) {
      $command .= ' TRAILING:' . $minQuality;
    }
    $command .= ' MINLEN:' . $minLength;
    $command .= ' 2> ' . $outputTrimLog;
    $command .= ' &&';
    $command .= ' grep \"^Input Read\" ' . $outputTrimLog . '| sed \'s/Input Read Pairs: \\([0-9]\\+\\).*Both Surviving: \\([0-9]\\+\\).*Forward Only Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\3/g\' | tr \'#\' \'\n\' > ' . $outputTrimStats;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

sub singleCommand {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rH_laneInfo = shift;
  my $outputDir   = shift;

  my $minQuality  = LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int');
  my $minLength   = LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int');
  my $adapterFile = LoadConfig::getParam($rH_cfg, 'trim', 'adapterFile', 1, 'filepath');
  my $headcrop    = LoadConfig::getParam($rH_cfg, 'trim', 'headcrop', 0, 'int');

  my $rawReadDir    = LoadConfig::getParam($rH_cfg, 'trim', 'rawReadDir', 1, 'dirpath');

  my $outputFastqName = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . $minQuality . 'l' . $minLength . '.single.fastq.gz';
  my $outputTrimLog = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.trim.out';
  my $outputTrimStats = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.trim.stats.csv';
  my $inputFastqName = $rawReadDir . '/' . $sampleName . '/run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '/' . $rH_laneInfo->{'read1File'};

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputFastqName], [$outputFastqName, $outputTrimLog, $outputTrimStats]);

  $ro_job->setOutputFileHash({SINGLE1_OUTPUT => $outputFastqName});

  if (!$ro_job->isUp2Date()) {
    my $command = "";
    $command .= LoadConfig::moduleLoad($rH_cfg, [
      ['trim', 'moduleVersion.java'],
      ['trim', 'moduleVersion.trimmomatic']
    ]);
    $command .= ' && java -XX:ParallelGCThreads=1 -Xmx2G -cp \$TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticSE';
    $command .= ' -threads ' . LoadConfig::getParam($rH_cfg, 'trim', 'nbThreads', 1, 'int');
    if ($rH_laneInfo->{'qualOffset'} eq "64") {
      $command .= ' -phred64';
    } else {
      $command .= ' -phred33';
    }
    $command .= ' ' . $inputFastqName;
    $command .= ' ' . $outputFastqName;
    if ($rH_laneInfo->{'qualOffset'} eq "64") {
      $command .= ' TOPHRED33';
    }
    if (defined($headcrop) && length($headcrop) > 0 && $headcrop > 0) {
      $command .= ' HEADCROP:' . $headcrop;
    }
    $command .= ' ILLUMINACLIP:' . $adapterFile . LoadConfig::getParam($rH_cfg, 'trim', 'clipSettings');
    if ($minQuality > 0) {
      $command .= ' TRAILING:' . $minQuality;
    }
    $command .= ' MINLEN:' . $minLength;
    $command .= ' 2> ' . $outputTrimLog;
    $command .= ' &&';
    $command .= ' grep \"^Input Read\" ' . $outputTrimLog . '| sed \'s/Input Reads: \\([0-9]\\+\\).*Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\2/g\' | tr \'#\' \'\n\' > ' . $outputTrimStats;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

1;
