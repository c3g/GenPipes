#!/usr/bin/env perl

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

  my $rO_job = new Job([$tumorBam, $normalBam], [$outputFile]);

  my $command;
  $rO_job->addModules($rH_cfg, [['countBins', 'moduleVersion.java'], ['countBins', 'moduleVersion.bvatools']]);
  $command .= "java -Djava.io.tmpdir=" . LoadConfig::getParam($rH_cfg, 'countBins', 'tmpDir') . " " . LoadConfig::getParam($rH_cfg, 'countBins', 'extraJavaFlags');
  $command .= " -Xmx1500M -jar \\\${BVATOOLS_JAR} bincounter";
  $command .= " \\\n  --norm " . $normType;
  $command .= " \\\n  --minMapQ " . LoadConfig::getParam($rH_cfg, 'countBins', 'minMapQ');
  $command .= " \\\n  --bam " . $tumorBam;
  if (defined($normalBam)) {
    $command .= " \\\n  --refbam " . $normalBam;
  }
  $command .= " \\\n  --window " . $window;
  $command .= " \\\n  > " . $outputFile;

  $rO_job->addCommand($command);

  return $rO_job;
}

sub deleteDuplicates {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $pair1           = shift;
  my $pair2           = shift;
  my $single          = shift;
  my $optOutputPrefix = shift;

  my $rO_job;
  if (defined($pair1) && defined($pair2)) {
    $rO_job = deletePairedDuplicates($rH_cfg, $sampleName, $pair1, $pair2, $optOutputPrefix);
  } elsif (defined($single)) {
    $rO_job = deleteSingleDuplicates($rH_cfg, $sampleName, $single, $optOutputPrefix);
  } else {
    die "Unknown runType. \n";
  }

  return $rO_job;
}

sub deletePairedDuplicates {
  my $rH_cfg       = shift;
  my $sampleName   = shift;
  my $pair1        = shift;
  my $pair2        = shift;
  my $outputPrefix = shift;
  my %retVal;

  my $command = "";
  my $outputFastqPair1Name = $outputPrefix . ".pair1.dup.fastq.gz";
  my $outputFastqPair2Name = $outputPrefix . ".pair2.dup.fastq.gz";

  my $rO_job = new Job([$pair1, $pair2], [$outputFastqPair1Name, $outputFastqPair2Name]);

  $rO_job->addModules($rH_cfg, [['default', 'moduleVersion.java'], ['duplicate', 'moduleVersion.bvatools']]);
  $command .= "java " . LoadConfig::getParam($rH_cfg, 'duplicate', 'extraJavaFlags') . " -Xmx" . LoadConfig::getParam($rH_cfg, 'duplicate', 'dupRam') . " -jar \\\$BVATOOLS_JAR filterdups" . " \\\n  --read1 " . $pair1 . " \\\n  --read2 " . $pair2;
  $command .= " \\\n  -k 20 -o 15";
  $command .= " \\\n  && mv " . $pair1 . ".dup.read1.gz " . $outputFastqPair1Name;
  $command .= " \\\n  && mv " . $pair2 . ".dup.read2.gz " . $outputFastqPair2Name;

  $rO_job->addCommand($command);

  return $rO_job;
}

sub deleteSingleDuplicates {
  my $rH_cfg       = shift;
  my $sampleName   = shift;
  my $single       = shift;
  my $outputPrefix = shift;
  my %retVal;

  my $outputFastqName = $outputPrefix . ".single.dup.fastq.gz";

  my $rO_job = new Job([$single], [$outputFastqName]);

  my $command;
  $rO_job->addModules($rH_cfg, [['default', 'moduleVersion.java'], ['duplicate', 'moduleVersion.bvatools']]);
  $command .= "java " . LoadConfig::getParam($rH_cfg, 'duplicate', 'extraJavaFlags') . " -Xmx" . LoadConfig::getParam($rH_cfg, 'duplicate', 'dupRam') . " -jar \\\$BVATOOLS_JAR filterdups" . " \\\n  --read1 " . $single;
  $command .= " \\\n  -k 20 -o 15";
  $command .= " \\\n  && mv " . $single . ".dup.read1.gz " . $outputFastqName;

  $rO_job->addCommand($command);

  return $rO_job;
}

sub resolveSampleBED {
  my $rH_cfg      = shift;
  my $readSet = shift;

  my $iniRegions = LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'coverageTargets', 0);

  if (!defined($iniRegions) || length($iniRegions) == 0) {
    return undef;
  } elsif ($iniRegions eq "auto") {
    if (@{$readSet->getBEDs()} == 0) {
      return undef;
    } else {
      return $readSet->getBEDs()->[0];
    }
  }

  return $iniRegions
}

sub depthOfCoverage {
  my $rH_cfg        = shift;
  my $inputBam      = shift;
  my $outputFile    = shift;
  my $coverageBED   = shift;
  my $refGenome     = shift;

  my $rO_job = new Job([$inputBam, $coverageBED], [$outputFile]);

  if (!defined($refGenome)) {
    $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta', 1, 'filepath');
  }
  
  my $maxDepth = LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'maxDepth', 0, 'array');

  my $command;
  $rO_job->addModules($rH_cfg, [['depthOfCoverage', 'moduleVersion.java'], ['depthOfCoverage', 'moduleVersion.bvatools']]);
  $command .= "java " . LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'extraJavaFlags') . " -Xmx" . LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'ram') . " -jar \\\${BVATOOLS_JAR}";
  $command .= " depthofcoverage";
  $command .= " \\\n  --simpleChrName";
  $command .= " \\\n  " . LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'extraFlags', 0);
  $command .= " \\\n  --ref " . $refGenome;
  if (defined($coverageBED) && length($coverageBED) > 0) {
    $command .= " \\\n  --intervals \\\'".$coverageBED . "'";
  }

  $command .= " \\\n  --bam " . $inputBam;
  $command .= " \\\n  > " . $outputFile;

  $rO_job->addCommand($command);

  return $rO_job;
}

sub qc {
  my $rH_cfg         = shift;
  my $read1          = shift;
  my $read2          = shift;
  my $type           = shift;
  my $regionName     = shift;
  my $outputDirectory= shift;
  
  my $rO_job = new Job([$read1, $read2], [$outputDirectory."/mpsQC_".$regionName."_stats.xml"]);
  
  my $nbThreads = LoadConfig::getParam($rH_cfg, 'generateQCGraphs','nbThreads');
  my $command;
  $rO_job->addModules($rH_cfg, [['generateQCGraphs', 'moduleVersion.java'], ['generateQCGraphs', 'moduleVersion.bvatools']]);
  $command .= "java " .LoadConfig::getParam($rH_cfg, 'generateQCGraphs', 'extraJavaFlags')." -Xmx".LoadConfig::getParam($rH_cfg, 'generateQCGraphs', 'maxRam')." -jar \\\${BVATOOLS_JAR}";
  $command .= ' readsqc --regionName \'' . $regionName . '\' --type ' . $type . ' --output \'' . $outputDirectory . '\' --read1 \'' . $read1 .'\'';
  if (defined($read2)) {
    $command .= ' --read2 \'' . $read2 . '\'';
  }
  if (defined($nbThreads)) {
    $command .= " --threads " . $nbThreads;
  }
  $rO_job->addCommand($command);

  return $rO_job;
}

1;
