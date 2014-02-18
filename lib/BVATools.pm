#!/usr/env/perl

=head1 NAME

I<BVATools>

=head1 SYNOPSIS

=head1 DESCRIPTION

B<SVTools> is a library to analyse BAMs for Structural Variants

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package BVATools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub countBins {
  my $rH_cfg     = shift;
  my $sampleName = shift;
  my $tumorBam   = shift;
  my $window     = shift;
  my $normType   = shift;
  my $outputFile = shift;
  my $normalBam  = shift; # can be undef for non paired run

  my $ro_job = new Job();
  if (defined($normalBam)) {
    $ro_job->testInputOutputs([$tumorBam, $normalBam], [$outputFile], $ro_job);
  } else {
    $ro_job->testInputOutputs([$tumorBam], [$outputFile], $ro_job);
  }

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['countBins', 'moduleVersion.java'], ['countBins', 'moduleVersion.bvatools']]) . ' && ';
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'countBins', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'countBins', 'extraJavaFlags');
    $command .= ' -Xmx1500M -jar \${BVATOOLS_JAR} bincounter';
    $command .= ' --norm ' . $normType;
    $command .= ' --minMapQ ' . LoadConfig::getParam($rH_cfg, 'countBins', 'minMapQ');
    $command .= ' --bam ' . $tumorBam;
    if (defined($normalBam)) {
      $command .= ' --refbam ' . $normalBam;
    }
    $command .= ' --window ' . $window;
    $command .= ' > ' . $outputFile;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub deleteDuplicates {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $pair1           = shift;
  my $pair2           = shift;
  my $single          = shift;
  my $optOutputPrefix = shift;

  my $ro_job;
  if (defined($pair1) && defined($pair2)) {
    $ro_job = deletePairedDuplicates($rH_cfg, $sampleName, $pair1, $pair2, $optOutputPrefix);
  } elsif (defined($single)) {
    $ro_job = deleteSingleDuplicates($rH_cfg, $sampleName, $single, $optOutputPrefix);
  } else {
    die "Unknown runType. \n";
  }

  return $ro_job;
}

sub deletePairedDuplicates {
  my $rH_cfg       = shift;
  my $sampleName   = shift;
  my $pair1        = shift;
  my $pair2        = shift;
  my $outputPrefix = shift;
  my %retVal;

  my $command = '';
  my $outputFastqPair1Name = $outputPrefix . '.pair1.dup.fastq.gz';
  my $outputFastqPair2Name = $outputPrefix . '.pair2.dup.fastq.gz';

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$pair1, $pair2], [$outputFastqPair1Name, $outputFastqPair2Name]);

  if (!$ro_job->isUp2Date()) {
    $command .= LoadConfig::moduleLoad($rH_cfg, [['default', 'moduleVersion.java'], ['duplicate', 'moduleVersion.bvatools']]) . ' &&';
    $command .= ' java ' . LoadConfig::getParam($rH_cfg, 'duplicate', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'duplicate', 'dupRam') . ' -jar \$BVATOOLS_JAR filterdups' . ' --read1 ' . $pair1 . ' --read2 ' . $pair2;
    $command .= ' -k 20 -o 15';
    $command .= ' && mv ' . $pair1 . '.dup.read1.gz ' . $outputFastqPair1Name;
    $command .= ' && mv ' . $pair2 . '.dup.read2.gz ' . $outputFastqPair2Name;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

sub deleteSingleDuplicates {
  my $rH_cfg       = shift;
  my $sampleName   = shift;
  my $single       = shift;
  my $outputPrefix = shift;
  my %retVal;

  my $outputFastqName = $outputPrefix . '.single.dup.fastq.gz';

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$single], [$outputFastqName]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['default', 'moduleVersion.java'], ['duplicate', 'moduleVersion.bvatools']]) . ' &&';
    $command .= ' java ' . LoadConfig::getParam($rH_cfg, 'duplicate', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'duplicate', 'dupRam') . ' -jar \$BVATOOLS_JAR filterdups' . ' --read1 ' . $single;
    $command .= ' -k 20 -o 15';
    $command .= ' && mv ' . $single . '.dup.read1.gz ' . $outputFastqName;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

sub resolveSampleBED {
  my $rH_cfg      = shift;
  my $rH_laneInfo = shift;

  my $iniRegions = LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'coverageTargets', 0);

  if (!defined($iniRegions) || length($iniRegions) == 0) {
    return undef;
  } elsif ($iniRegions eq 'auto') {
    if (@{$rH_laneInfo->{'bedFiles'}} == 0) {
      return undef;
    } else {
      return $rH_laneInfo->{'bedFiles'}->[0];
    }
  }

  return $iniRegions
}

sub depthOfCoverage {
  my $rH_cfg      = shift;
  my $inputBam    = shift;
  my $outputFile  = shift;
  my $coverageBED = shift;

  my $ro_job = new Job();
  if (defined($coverageBED)) {
    $ro_job->testInputOutputs([$inputBam, $coverageBED], [$outputFile]);
  } else {
    $ro_job->testInputOutputs([$inputBam], [$outputFile]);
  }

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta', 1, 'filepath');
  my $rA_thresholds = LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'percentThresholds', 0, 'array');

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['depthOfCoverage', 'moduleVersion.java'], ['depthOfCoverage', 'moduleVersion.bvatools']]) . ' &&';
    $command .= ' java ' . LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'ram') . ' -jar \${BVATOOLS_JAR}';
    $command .= ' depthofcoverage';
    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'extraFlags', 0);
    $command .= ' --ref ' . $refGenome;
    if (defined($coverageBED) && length($coverageBED) > 0) {
      $command .= ' --intervals ' . $coverageBED;
    }

    if (defined($rA_thresholds) and $rA_thresholds ne "") {
      for my $threshold (@{$rA_thresholds}) {
        $command .= ' --summaryCoverageThresholds ' . $threshold;
      }
    }
    $command .= ' --bam ' . $inputBam;
    $command .= ' > ' . $outputFile;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

1;
