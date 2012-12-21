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
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub realign {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $seqName         = shift;
  my $processUnmapped = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $sortedBAM = $sampleName.'/'.$sampleName.'.sorted.bam';
  my $intervalOutput = $sampleName.'/realign/'.$seqName.'.intervals';
  my $realignOutput = $sampleName.'/realign/'.$seqName.'.bam';
  
  my $command;
  $command .= 'module load mugqic/GenomeAnalysisTKLite/2.1-13 ;';
  $command .= ' java '.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'realignRam').'  -jar \${GATK_JAR}';
  $command .= ' -T RealignerTargetCreator';
  $command .= ' -R '.$refGenome;
  $command .= ' -o '.$intervalOutput;
  $command .= ' -I '.$sortedBAM;
  $command .= ' -L '.$seqName;
  $command .= ' ; ';
  $command .= ' java '.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'indelRealigner', 'realignRam').' -jar \${GATK_JAR}';
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
  my $rH_cfg          = shift;
  my $sampleName      = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $sortedBAM = $sampleName.'/'.$sampleName.'.sorted.dup.bam';
  my $output = $sampleName.'.sorted.dup.coverage';
  
  my $command;
  $command .= 'module load mugqic/GenomeAnalysisTKLite/2.1-13 ;';
  $command .= ' java '.LoadConfig::getParam($rH_cfg, 'genomeCoverage', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'genomeCoverage', 'genomeCoverageRam').'  -jar \${GATK_JAR}';
  $command .= ' -T DepthOfCoverage --omitDepthOutputAtEachBase --logging_level ERROR';
  $command .= ' --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 75 --summaryCoverageThreshold 100';
  $command .= ' --start 1 --stop 1000 --nBins 999 -dt NONE';
  $command .= ' -R '.$refGenome;
  $command .= ' -o '.$output;
  $command .= ' -I '.$sortedBAM;
  
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
  $command .= 'module load mugqic/GenomeAnalysisTKLite/2.1-13 ;';
  $command .= ' java '.LoadConfig::getParam($rH_cfg, 'targetCoverage', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'targetCoverage', 'coverageRam').'  -jar \${GATK_JAR}';
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
