#!/usr/env/perl

=head1 NAME

I<Picard>

=head1 SYNOPSIS

Picard->merge()

=head1 DESCRIPTION

B<Picard> is a library to manipulate and compute stats off of BAMs

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
sub mergeFiles {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $rA_inputFiles = shift;
  my $outputBAM    = shift;

  my $latestBam;
  my $bamInputs;
  my $countInputs;
  for my $file (@$rA_inputFiles) {
    if(!defined($latestBam)) {
      $latestBam = -M $file;
    }
    else {
      my $modDate = -M $file;
      if(!defined($modDate)) {
        warn "Input merge file doesn't exist: $file\n";
      }
      if($modDate < $latestBam) {
        $latestBam = $modDate;
      }
    }
    $bamInputs .= 'INPUT='.$file.' ';
  }

  my $command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestBam) || !defined(-M $outputBAM) || $latestBam < -M $outputBAM) {
  $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'mergeFiles','moduleVersion.java') .' ' .LoadConfig::getParam($rH_cfg, 'mergeFiles', 'moduleVersion.picard').' &&'; 
    $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'mergeFiles', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'mergeFiles', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'mergeFiles', 'mergeRam').' -jar \${PICARD_HOME}/MergeSamFiles.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true';
    $command .= ' TMP_DIR='.LoadConfig::getParam($rH_cfg, 'mergeFiles', 'tmpDir');
    $command .= ' '.$bamInputs;
    $command .= ' OUTPUT='.$outputBAM;
    $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'mergeFiles', 'mergeRecInRam');
  }
  return $command;
}

sub fixmate {
  my $rH_cfg     = shift;
  my $inputBAM   = shift;
  my $outputBAM  = shift;

  my $command;
  $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'fixmate','moduleVersion.java') .' ' .LoadConfig::getParam($rH_cfg, 'fixmate', 'moduleVersion.picard').' &&'; 
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'fixmate', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'fixmate', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'fixmate', 'fixmateRam').' -jar \${PICARD_HOME}/FixMateInformation.jar';
  $command .= ' VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate';
  $command .= ' TMP_DIR='.LoadConfig::getParam($rH_cfg, 'fixmate', 'tmpDir');
  $command .= ' INPUT='.$inputBAM;
  $command .= ' OUTPUT='.$outputBAM;
  $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'fixmate', 'fixmateRecInRam');
  return $command;
}

sub markDup {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBAM      = shift;
  my $outputBAM     = shift;
  my $outputMetrics = shift;

#   if (!(defined $inputBAM)) {
#     $inputBAM = $sampleName.'/'.$sampleName.'.matefixed.sorted.bam';
#   }
#   if (!(defined $outputBAM)) {
#     $outputBAM = $sampleName.'/'.$sampleName.'.sorted.dup.bam';
#   }
#   if (!(defined $outputMetrics)) {
#     $outputMetrics = $sampleName.'/'.$sampleName.'.sorted.dup.metrics';
#   }

  my $command;
  $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'markDup','moduleVersion.java') .' ' .LoadConfig::getParam($rH_cfg, 'markDup', 'moduleVersion.picard').' &&'; 
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'markDup', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'markDupRam', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'markDup', 'markDupRam').' -jar \${PICARD_HOME}/MarkDuplicates.jar';
  $command .= ' REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true';
  $command .= ' TMP_DIR='.LoadConfig::getParam($rH_cfg, 'markDup', 'tmpDir');
  $command .= ' INPUT='.$inputBAM;
  $command .= ' OUTPUT='.$outputBAM;
  $command .= ' METRICS_FILE='.$outputMetrics;
  $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'markDup', 'markDupRecInRam');

  return $command;
}

sub collectMetrics {
  my $rH_cfg        = shift;
  my $inputBAM      = shift;
  my $outputMetrics = shift;

  my $command;
  $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'collectMetrics','moduleVersion.java') .' ' .LoadConfig::getParam($rH_cfg, 'collectMetrics', 'moduleVersion.picard').' &&'; 
  $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'collectMetricsRam').' -jar \${PICARD_HOME}/CollectMultipleMetrics.jar';
  $command .= ' PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics  VALIDATION_STRINGENCY=SILENT';
  $command .= ' TMP_DIR='.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'tmpDir');
  $command .= ' REFERENCE_SEQUENCE='.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'referenceFasta');
  $command .= ' INPUT='.$inputBAM;
  $command .= ' OUTPUT='.$outputMetrics;
  $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'collectMetricsRecInRam');

  return $command;
}

# Sort BAM/SAM files
sub sortSam {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBAM    = shift;
  my $outputBAM     = shift;
  my $order         = shift;

  my $latestBam = -M $inputBAM;

  my $command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestBam) || !defined(-M $outputBAM) || $latestBam < -M $outputBAM) {
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'sortSam','moduleVersion.java') .' ' .LoadConfig::getParam($rH_cfg, 'sortSam', 'moduleVersion.picard') .' &&';
    $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'sortSam', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'sortSam', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'sortSam', 'sortRam').' -jar \${PICARD_HOME}/SortSam.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true';
    $command .= ' TMP_DIR='.LoadConfig::getParam($rH_cfg, 'sortSam', 'tmpDir');
    $command .= ' INPUT='.$inputBAM;
    $command .= ' OUTPUT='.$outputBAM;
    $command .= ' SORT_ORDER='.$order;
    $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'sortSam', 'sortRecInRam');
  }
  return $command;
}


# reorder BAM/SAM files based on reference/dictionary
sub reorderSam {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBAM      = shift;
  my $outputBAM     = shift;


  my $latestBam = -M $inputBAM;

  my $command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestBam) || !defined(-M $outputBAM) || $latestBam < -M $outputBAM) {
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'reorderSam','moduleVersion.java') .' ' .LoadConfig::getParam($rH_cfg, 'reorderSam', 'moduleVersion.picard') .' &&';
    $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'reorderSam', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'reorderSam', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'reorderSam', 'reorderRam').' -jar \${PICARD_HOME}/ReorderSam.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true';
    $command .= ' TMP_DIR='.LoadConfig::getParam($rH_cfg, 'reorderSam', 'tmpDir');
    $command .= ' INPUT='.$inputBAM;
    $command .= ' OUTPUT='.$outputBAM;
    $command .= ' REFERENCE='.LoadConfig::getParam($rH_cfg, 'reorderSam', 'referenceFasta');
    $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'reorderSam', 'reorderRecInRam');
  }
  return $command;
}

1;
