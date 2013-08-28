#!/usr/env/perl

=head1 NAME

I<GATK>

=head1 SYNOPSIS

GATK->trim()

=head1 DESCRIPTION

B<GATK> is a library to manipulate and compute stats off of BAMs

Input = file_name

Output = array

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package GATK;

# Strict Pragmas
#--------------------------
use PipelineUtils;
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub recalibration {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $sortedBAM       = shift;
  my $outputPrefix    = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $knownSites = LoadConfig::getParam($rH_cfg, 'recalibration', 'knownSites');
  my $recalOutput = $outputPrefix.'.recalibration_report.grp';
  my $bamOutput = $outputPrefix.'.recal.bam';

  my $command;
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'recalibration', 'moduleVersion.java').' '.LoadConfig::getParam($rH_cfg, 'recalibration', 'moduleVersion.gatk').' ;';
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'recalibration', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'recalibration', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'recalibration', 'recalRam').'  -jar \${GATK_JAR}';
  $command .= ' -T BaseRecalibrator';
  $command .= ' -nct '.LoadConfig::getParam($rH_cfg, 'recalibration', 'threads');
  $command .= ' -R '.$refGenome;
  $command .= ' -knownSites '.$knownSites;
  $command .= ' -o '.$recalOutput;
  $command .= ' -I '.$sortedBAM;
  $command .= ' ; ';
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'recalibration', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'recalibration', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'recalibration', 'recalRam').' -jar \${GATK_JAR}';
  $command .= ' -T PrintReads';
  $command .= ' -nct '.LoadConfig::getParam($rH_cfg, 'recalibration', 'threads');
  $command .= ' -R '.$refGenome;
  $command .= ' -BQSR '.$recalOutput;
  $command .= ' -o '.$bamOutput;
  $command .= ' -I '.$sortedBAM;

  return $command;
}

sub realign {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $sortedBAM       = shift;
  my $seqName         = shift;
  my $outputPrefix    = shift;
  my $processUnmapped = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $intervalOutput = $outputPrefix.'/'.$sampleName.'/realign/'.$seqName.'.intervals';
  my $realignOutput = $outputPrefix.'/'.$sampleName.'/realign/'.$seqName.'.bam';
  
  my $command;
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'moduleVersion.java').' '.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'moduleVersion.gatk').' ;';
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'realignRam').'  -jar \${GATK_JAR}';
  $command .= ' -T RealignerTargetCreator';
  $command .= ' -R '.$refGenome;
  $command .= ' -o '.$intervalOutput;
  $command .= ' -I '.$sortedBAM;
  $command .= ' -L '.$seqName;
  $command .= ' ; ';
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'realignRam').' -jar \${GATK_JAR}';
  $command .= ' -T IndelRealigner';
  $command .= ' -R '.$refGenome;
  $command .= ' -targetIntervals '.$intervalOutput;
  $command .= ' -o '.$realignOutput;
  $command .= ' -I '.$sortedBAM;
  $command .= ' -L '.$seqName;
  $command .= ' --maxReadsInMemory '.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'realignReadsInRam');
  if($processUnmapped == 1) {
    $command .= ' -L unmapped';
  }
  
  return $command;
}

sub genomeCoverage {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBam      = shift;
  my $outputPrefix  = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $rA_thresholds = LoadConfig::getParam($rH_cfg, 'genomeCoverage', 'percentThresholds');
  
  my $command;
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'genomeCoverage', 'moduleVersion.java').' '.LoadConfig::getParam($rH_cfg, 'genomeCoverage', 'moduleVersion.gatk').' ;';
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'genomeCoverage', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'genomeCoverage', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'genomeCoverage', 'genomeCoverageRam').'  -jar \${GATK_JAR}';
  $command .= ' -T DepthOfCoverage --omitDepthOutputAtEachBase --logging_level ERROR';
  my $highestThreshold = 0;
  for my $threshold (@$rA_thresholds) {
    $command .= ' --summaryCoverageThreshold '.$threshold;
    if($highestThreshold > $threshold) {
      die "Tresholds must be ascending: ".join(',', @$rA_thresholds);
    }
    $highestThreshold = $threshold;
  }
  $command .= ' --start 1 --stop '.$highestThreshold.' --nBins '.($highestThreshold-1).' -dt NONE';
  $command .= ' -R '.$refGenome;
  $command .= ' -o '.$outputPrefix;
  $command .= ' -I '.$inputBam;
  
  return $command;
}

sub targetCoverage {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBam      = shift;
  my $outputPrefix  = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $targets = LoadConfig::getParam($rH_cfg, 'targetCoverage', 'coverageTargets');
  my $rA_thresholds = LoadConfig::getParam($rH_cfg, 'targetCoverage', 'percentThresholds');

  my $command = "";
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'targetCoverage', 'moduleVersion.java').' '.LoadConfig::getParam($rH_cfg, 'targetCoverage', 'moduleVersion.gatk').' ;';
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'targetCoverage', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'targetCoverage', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'targetCoverage', 'coverageRam').'  -jar \${GATK_JAR}';
  $command .= ' -T DepthOfCoverage --omitDepthOutputAtEachBase --logging_level ERROR';
  my $highestThreshold = 0;
  for my $threshold (@$rA_thresholds) {
    $command .= ' --summaryCoverageThreshold '.$threshold;
    if($highestThreshold > $threshold) {
      die "Tresholds must be ascending: ".join(',', @$rA_thresholds);
    }
    $highestThreshold = $threshold;
  }
  $command .= ' --start 1 --stop '.$highestThreshold.' --nBins '.($highestThreshold-1).' -dt NONE';
  $command .= ' -R '.$refGenome;
  $command .= ' -o '.$outputPrefix;
  $command .= ' -I '.$inputBam;
  $command .= ' -L '.$targets;
  
  return $command;
}

1;
