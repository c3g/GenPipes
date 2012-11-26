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

package Picard;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub merge {
  my $rH_cfg            = shift;
  my $sampleName        = shift;
  my $rAoH_sampleLanes  = shift;

  my $latestBam;
  my $bamInputs;
  my $countInputs;
  my $outputBAM = $sampleName.'/'.$sampleName.'.sorted.bam';
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $directory = $sampleName."/run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'}."/";
    my $sortedLaneBamFile = $directory.$rH_laneInfo->{'name'}.".sorted.bam";
    my $laneStatsFile = $directory.$rH_laneInfo->{'name'}.".counts";
    my $runName = $sampleName."_run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'};


    if(!defined($latestBam)) {
      $latestBam = -M $sortedLaneBamFile;
    }
    else {
      my $modDate = -M $sortedLaneBamFile;
      if($modDate < $latestBam) {
        $latestBam = $modDate;
      }
    }
    $bamInputs .= 'INPUT='.$sortedLaneBamFile.' ';
    $countInputs .= $laneStatsFile.' ';
  }

  my $command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestBam) || !defined(-M $outputBAM) || $latestBam < -M $outputBAM) {
    $command .= 'module load mugqic/picard/1.77 ;';
    $command .= ' java -Xmx'.LoadConfig::getParam($rH_cfg, 'mergeLanes', 'mergeRam').' -jar \${PICARD_HOME}/MergeSamFiles.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true';
    $command .= ' '.$bamInputs;
    $command .= ' OUTPUT='.$outputBAM;
    $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'mergeLanes', 'mergeRecInRam');
    $command .= ' ; cat '.$countInputs.' > $sampleName/$sampleName.runLane.counts';
  }
  return $command;
}

sub mergeRealigned {
  my $rH_cfg            = shift;
  my $sampleName        = shift;
  my $rA_seqNames  = shift;

  my $latestBam;
  my $bamInputs;
  my $countInputs;
  my $outputBAM = $sampleName.'/'.$sampleName.'.realigned.qsorted.bam';
  for my $seqName (@$rA_seqNames) {
    my $directory = $sampleName."/realign/";
    my $realignedBamFile = $directory.$seqName.'.bam';

    $bamInputs .= 'INPUT='.$realignedBamFile.' ';
  }

  my $command;
  $command .= 'module load mugqic/picard/1.77 ;';
  $command .= ' java -Xmx'.LoadConfig::getParam($rH_cfg, 'mergeRealigned', 'mergeRam').' -jar \${PICARD_HOME}/MergeSamFiles.jar';
  $command .= ' VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=false CREATE_INDEX=true SORT_ORDER=queryname';
  $command .= ' '.$bamInputs;
  $command .= ' OUTPUT='.$outputBAM;
  $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'mergeRealigned', 'mergeRecInRam');
  return $command;
}

sub fixmate {
  my $rH_cfg     = shift;
  my $sampleName = shift;

  my $inputBAM = $sampleName.'/'.$sampleName.'.realigned.qsorted.bam';
  my $outputBAM = $sampleName.'/'.$sampleName.'.matefixed.sorted.bam';

  my $command;
  $command .= 'module load mugqic/picard/1.77 ;';
  $command .= ' java -Xmx'.LoadConfig::getParam($rH_cfg, 'fixmate', 'fixmateRam').' -jar \${PICARD_HOME}/FixMateInformation.jar';
  $command .= ' VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate';
  $command .= ' INPUT='.$inputBAM;
  $command .= ' OUTPUT='.$outputBAM;
  $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'fixmate', 'fixmateRecInRam');
  return $command;
}

sub markDup {
  my $rH_cfg     = shift;
  my $sampleName = shift;

  my $inputBAM = $sampleName.'/'.$sampleName.'.matefixed.sorted.bam';
  my $outputBAM = $sampleName.'/'.$sampleName.'.sorted.dup.bam';
  my $outputMetrics = $sampleName.'/'.$sampleName.'.sorted.dup.metrics';
  
  my $command;
  $command .= 'module load mugqic/picard/1.77 ;';
  $command .= ' java -Xmx'.LoadConfig::getParam($rH_cfg, 'markDup', 'markDupRam').' -jar \${PICARD_HOME}/MarkDuplicates.jar';
  $command .= ' REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true';
  $command .= ' INPUT='.$inputBAM;
  $command .= ' OUTPUT='.$outputBAM;
  $command .= ' METRICS_FILE='.$outputMetrics;
  $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'markDup', 'markDupRecInRam');

  return $command;
}

sub collectMetrics {
  my $rH_cfg     = shift;
  my $sampleName = shift;

  my $inputBAM = $sampleName.'/'.$sampleName.'.sorted.dup.bam';
  my $outputMetrics = $sampleName.'/'.$sampleName.'.sorted.dup.all.metrics';
  
  my $command;
  $command .= 'module load mugqic/picard/1.77 ;';
  $command .= ' java -Xmx'.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'collectMetricsRam').' -jar \${PICARD_HOME}/CollectMultipleMetrics.jar';
  $command .= ' PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics  VALIDATION_STRINGENCY=SILENT';
  $command .= ' REFERENCE_SEQUENCE='.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'referenceFasta');
  $command .= ' INPUT='.$inputBAM;
  $command .= ' OUTPUT='.$outputMetrics;
  $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'collectMetricsRecInRam');

  return $command;
}
1;
