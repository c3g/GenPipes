#!/usr/bin/env perl

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
use File::Basename;
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
#  my $sampleName  = shift;
#  my $rH_laneInfo = shift;
  my $input1 = shift;   # Contains BAM or FASTQ pair1 or FASTQ single
  my $input2 = shift;   # Contains FASTQ pair2 if any
  my $outputDir  = shift;    # FASTQ output file(s) prefix including directories if any
  my $qualOffset = shift;

  my $ro_job = new Job();

#  my $rawReadDir = LoadConfig::getParam($rH_cfg, 'trim', 'rawReadDir', 1, 'dirpath');
#  my $inputDir = "$rawReadDir/$sampleName/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";

#  my $inputFastqPair1Name = $inputDir . $rH_laneInfo->{'read1File'};
#  my $inputFastqPair2Name = undef;

#  my $skipTrimming = LoadConfig::getParam($rH_cfg, 'trim', 'skip', 0);
#
#  if ($rH_laneInfo->{'runType'} eq "SINGLE_END") {
#    if ($skipTrimming) {
#      $ro_job->setOutputFileHash({SINGLE1_OUTPUT => $inputFastqPair1Name});
#      $ro_job->setUp2Date(1);
#    } else {
#      $ro_job = singleCommand($rH_cfg, $sampleName, $rH_laneInfo, $outputDir);
#    }
#  } elsif ($rH_laneInfo->{'runType'} eq "PAIRED_END") {
#    if ($skipTrimming) {
#      $inputFastqPair2Name = $inputDir . $rH_laneInfo->{'read2File'};
#      $ro_job->setOutputFileHash({SINGLE1_OUTPUT => $inputFastqPair1Name, PAIR1_OUTPUT => $inputFastqPair1Name, PAIR2_OUTPUT => $inputFastqPair2Name});
#      $ro_job->setUp2Date(1);
#    } else {
#      $ro_job = pairCommand($rH_cfg, $sampleName, $rH_laneInfo, $outputDir);
#    }
#  } else {
#    die "Error in Trimmomatic::trim: unknown runType: " . $rH_laneInfo->{'runType'} . "!";
#  }

  my $minQuality  = LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int');
  my $minLength   = LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int');

  my $outputPrefix = $outputDir . "/" . basename($input1);
  $outputPrefix =~ s/\.(single|pair[12]).fastq(.gz)?$//;
  my $outputFilePrefix = $outputPrefix . ".t" . $minQuality . "l" . $minLength;

  my $inputs;
  my $outputs;

  if ($input2) {    # Paired end reads
    $inputs = [$input1, $input2];
    $outputs = [
      $outputFilePrefix . '.pair1.fastq.gz',
      $outputFilePrefix . '.single1.fastq.gz',
      $outputFilePrefix . '.pair2.fastq.gz',
      $outputFilePrefix . '.single2.fastq.gz'
    ];
  } else {    # Single end reads
    $inputs = [$input1];
    $outputs = [$outputFilePrefix . '.single.fastq.gz'];
  }

  my $outputTrimLog = "$outputPrefix.trim.out";
  my $outputTrimStats = "$outputPrefix.trim.stats.csv";

  push(@$outputs, $outputTrimLog, $outputTrimStats);

  $ro_job->testInputOutputs($inputs, $outputs);

  if (!$ro_job->isUp2Date()) {
    my $command = "mkdir -p " . dirname($outputPrefix) . " && \\\n";

    $command .= LoadConfig::moduleLoad($rH_cfg, [
      ['trim', 'moduleVersion.java'],
      ['trim', 'moduleVersion.trimmomatic']
    ]) . " && \\\n";

    $command .= "java -XX:ParallelGCThreads=1 -Xmx2G -cp \\\$TRIMMOMATIC_JAR " . ($input2 ? "PE" : "SE") . " \\\n";
    $command .= "  -threads " . LoadConfig::getParam($rH_cfg, 'trim', 'nbThreads', 1, 'int') . " \\\n";
    $command .= "  -phred" . ($qualOffset eq "64" ? $qualOffset : "33") . " \\\n";
    $command .= "  " . join(" \\\n  ", @$inputs) . " \\\n";
    $command .= "  " . join(" \\\n  ", @$outputs) . " \\\n";
    $command .= "  ILLUMINACLIP:" . LoadConfig::getParam($rH_cfg, 'trim', 'adapterFile', 1, 'filepath') . LoadConfig::getParam($rH_cfg, 'trim', 'clipSettings') . " \\\n";
    if ($minQuality > 0) {
      $command .= "  TRAILING:$minQuality \\\n";
    }
    my $headcrop    = LoadConfig::getParam($rH_cfg, 'trim', 'headcrop', 0, 'int');
    if ($headcrop and $headcrop > 0) {
      $command .= "  HEADCROP:$headcrop \\\n";
    }
    $command .= "  MINLEN:$minLength \\\n";
    $command .= $qualOffset eq "64" ? "  TOPHRED33 \\\n" : "";
    $command .= "  2> $outputTrimLog && \\\n";

    # Compute trim stats
    if ($input2) {
      $command .= ' grep \"^Input Read\" ' . $outputTrimLog . '| sed \'s/Input Read Pairs: \\([0-9]\\+\\).*Both Surviving: \\([0-9]\\+\\).*Forward Only Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\3/g\' | tr \'#\' \'\n\' > ' . $outputTrimStats . " \\\n";
    } else {
      $command .= ' grep \"^Input Read\" ' . $outputTrimLog . '| sed \'s/Input Reads: \\([0-9]\\+\\).*Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\2/g\' | tr \'#\' \'\n\' > ' . $outputTrimStats;
    }

    $ro_job->addCommand($command);
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
  my $headcrop    = LoadConfig::getParam($rH_cfg, 'trim', 'headcrop', 0, 'int');
  my $rawReadDir  = LoadConfig::getParam($rH_cfg, 'trim', 'rawReadDir', 1, 'dirpath');

  my $outputPrefix = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . $minQuality . 'l' . $minLength;
  my $outputFastqPair1Name   = $outputPrefix . '.pair1.fastq.gz';
  my $outputFastqPair2Name   = $outputPrefix . '.pair2.fastq.gz';
  my $outputFastqSingle1Name = $outputPrefix . '.single1.fastq.gz';
  my $outputFastqSingle2Name = $outputPrefix . '.single2.fastq.gz';
  my $outputTrimLog = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.trim.out';
  my $outputTrimStats = $outputDir . '/' . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.trim.stats.csv';

  my $inputDir = $rawReadDir . '/' . $sampleName . '/run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '/';

  my $inputBam;
  my $inputFastqPair1Name;
  my $inputFastqPair2Name;

  my $rawReadFormat = LoadConfig::getParam($rH_cfg, 'default', 'rawReadFormat', 0);
  my $rA_inputs;

  if (defined($rawReadFormat) and $rawReadFormat eq "fastq") {
    $inputFastqPair1Name = $inputDir . $rH_laneInfo->{'read1File'};
    $inputFastqPair2Name = $inputDir . $rH_laneInfo->{'read2File'};

    $rA_inputs = [$inputFastqPair1Name, $inputFastqPair2Name];
  } else {    # BAM format by default
    $inputBam = $inputDir . $rH_laneInfo->{'read1File'};
    $inputFastqPair1Name = $inputBam;
    $inputFastqPair1Name =~ s/bam$/pair1.fastq.gz/;
    $inputFastqPair2Name = $inputBam;
    $inputFastqPair2Name =~ s/bam$/pair2.fastq.gz/;
    $rA_inputs = [$inputBam];
  }

  my $ro_job = new Job();
  $ro_job->testInputOutputs($rA_inputs, [$outputFastqPair1Name, $outputFastqPair2Name, $outputFastqSingle1Name, $outputFastqSingle2Name, $outputTrimLog, $outputTrimStats]);

  $ro_job->setOutputFileHash({PAIR1_OUTPUT => $outputFastqPair1Name, PAIR2_OUTPUT => $outputFastqPair2Name, SINGLE1_OUTPUT => $outputFastqSingle1Name, SINGLE2_OUTPUT => $outputFastqSingle2Name});

  if (!$ro_job->isUp2Date()) {
    my $command = "mkdir -p $outputDir && \\\n";

    $command .= LoadConfig::moduleLoad($rH_cfg, [
      ['trim', 'moduleVersion.java'],
      ['trim', 'moduleVersion.trimmomatic']
    ]) . " && \\\n";

    unless (defined($rawReadFormat) and $rawReadFormat eq "fastq") {   # Assume BAM format by default
      # Convert BAM file to paired FASTQ files
      $command .= LoadConfig::moduleLoad($rH_cfg, [
        ['trim', 'moduleVersion.picard']
      ]) . " && \\\n";
      $command .= "java -Djava.io.tmpdir=" . LoadConfig::getParam($rH_cfg, 'trim', 'tmpDir') . " " . LoadConfig::getParam($rH_cfg, 'trim', 'extraJavaFlags') . " -jar \\\$PICARD_HOME/SamToFastq.jar \\\n";
      $command .= "  INPUT=$inputBam \\\n";
      $command .= "  FASTQ=$inputFastqPair1Name \\\n";
      $command .= "  SECOND_END_FASTQ=$inputFastqPair2Name && \\\n";
    }

    $command .= "java -XX:ParallelGCThreads=1 -Xmx2G -cp \\\$TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE \\\n";
    $command .= "  -threads " . LoadConfig::getParam($rH_cfg, 'trim', 'nbThreads', 1, 'int') . " \\\n";
    if ($rH_laneInfo->{'qualOffset'} eq "64") {
      $command .= "  -phred64 \\\n";
    } else {
      $command .= "  -phred33 \\\n";
    }
    $command .= "  $inputFastqPair1Name \\\n";
    $command .= "  $inputFastqPair2Name \\\n";
    $command .= "  $outputFastqPair1Name $outputFastqSingle1Name \\\n";
    $command .= "  $outputFastqPair2Name $outputFastqSingle2Name \\\n";
    if ($rH_laneInfo->{'qualOffset'} eq "64") {
      $command .= "  TOPHRED33 \\\n";
    }
    if (defined($headcrop) && length($headcrop) > 0 && $headcrop > 0) {
      $command .= "  HEADCROP:$headcrop \\\n";
    }
    $command .= "  ILLUMINACLIP:$adapterFile" . LoadConfig::getParam($rH_cfg, 'trim', 'clipSettings') . " \\\n";
    if ($minQuality > 0) {
      $command .= "  TRAILING:$minQuality \\\n";
    }
    $command .= "  MINLEN:$minLength \\\n";
    $command .= "  2> $outputTrimLog && \\\n";
    $command .= ' grep \"^Input Read\" ' . $outputTrimLog . '| sed \'s/Input Read Pairs: \\([0-9]\\+\\).*Both Surviving: \\([0-9]\\+\\).*Forward Only Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\3/g\' | tr \'#\' \'\n\' > ' . $outputTrimStats . " \\\n";

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
  # Input can be in FASTQ or BAM format
  $ro_job->testInputOutputs([$inputFastqName], [$outputFastqName, $outputTrimLog, $outputTrimStats]);

  $ro_job->setOutputFileHash({SINGLE1_OUTPUT => $outputFastqName});

  if (!$ro_job->isUp2Date()) {
    my $command = "";

    $command .= LoadConfig::moduleLoad($rH_cfg, [
      ['trim', 'moduleVersion.java'],
      ['trim', 'moduleVersion.trimmomatic']
    ]);

    my $rawReadFormat = LoadConfig::getParam($rH_cfg, 'default', 'rawReadFormat', 0);
    unless (defined($rawReadFormat) and $rawReadFormat eq "fastq") {   # Assume BAM format by default
      # Convert BAM file to single FASTQ files
      $command .= ' && ' . LoadConfig::moduleLoad($rH_cfg, [
        ['trim', 'moduleVersion.picard']
      ]);
      $command .= ' && java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'trim', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'trim', 'extraJavaFlags') . ' -jar \${PICARD_HOME}/SamToFastq.jar';
      $command .= ' INPUT=' . $inputFastqName;
      $inputFastqName =~ s/bam$/single.fastq.gz/;
      $command .= ' FASTQ=' . $inputFastqName;
    }

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
