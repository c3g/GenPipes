#!/usr/bin/env perl

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

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin";

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------
sub mergeFiles {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $rA_inputFiles = shift;
  my $outputBAM     = shift;

  my $bamInputs;
  my $countInputs;
  for my $file (@$rA_inputFiles) {
    $bamInputs .= 'INPUT=' . $file . ' ';
  }

  my $rO_job = new Job();
  $rO_job->testInputOutputs($rA_inputFiles, [$outputBAM]);

  if (!$rO_job->isUp2Date()) {
    my $command;
    $rO_job->addModules($rH_cfg, [['mergeFiles', 'moduleVersion.java'], ['mergeFiles', 'moduleVersion.picard']]); 
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'mergeFiles', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'mergeFiles', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'mergeFiles', 'mergeRam') . ' -jar \${PICARD_HOME}/MergeSamFiles.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true';
    $command .= ' TMP_DIR=' . LoadConfig::getParam($rH_cfg, 'mergeFiles', 'tmpDir');
    $command .= ' ' . $bamInputs;
    $command .= ' OUTPUT=' . $outputBAM;
    $command .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam($rH_cfg, 'mergeFiles', 'mergeRecInRam', 1, 'int');

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub fixmate {
  my $rH_cfg     = shift;
  my $inputBAM   = shift;
  my $outputBAM  = shift;

  my $rO_job = new Job();
  $rO_job->testInputOutputs([$inputBAM], [$outputBAM]);

  if (!$rO_job->isUp2Date()) {
    my $command;
    $rO_job->addModules($rH_cfg, [['fixmate', 'moduleVersion.java'], ['fixmate', 'moduleVersion.picard']]); 
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'fixmate', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'fixmate', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'fixmate', 'fixmateRam') . ' -jar \${PICARD_HOME}/FixMateInformation.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate';
    $command .= ' TMP_DIR=' . LoadConfig::getParam($rH_cfg, 'fixmate', 'tmpDir');
    $command .= ' INPUT=' . $inputBAM;
    $command .= ' OUTPUT=' . $outputBAM;
    $command .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam($rH_cfg, 'fixmate', 'fixmateRecInRam', 1, 'int');

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub markDup {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBAM      = shift;
  my $outputBAM     = shift;
  my $outputMetrics = shift;

#   if (!(defined $inputBAM)) {
#     $inputBAM = $sampleName . '/' . $sampleName . '.matefixed.sorted.bam';
#   }
#   if (!(defined $outputBAM)) {
#     $outputBAM = $sampleName . '/' . $sampleName . '.sorted.dup.bam';
#   }
#   if (!(defined $outputMetrics)) {
#     $outputMetrics = $sampleName . '/' . $sampleName . '.sorted.dup.metrics';
#   }
  my $rO_job = new Job();
  $rO_job->testInputOutputs([$inputBAM], [$outputBAM, $outputMetrics]);

  if (!$rO_job->isUp2Date()) {
    my $command;
    $rO_job->addModules($rH_cfg, [['markDup', 'moduleVersion.java'], ['markDup', 'moduleVersion.picard']]); 
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'markDup', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'markDup', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'markDup', 'markDupRam') . ' -jar \${PICARD_HOME}/MarkDuplicates.jar';
    $command .= ' REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true';
    $command .= ' TMP_DIR=' . LoadConfig::getParam($rH_cfg, 'markDup', 'tmpDir');
    $command .= ' INPUT=' . $inputBAM;
    $command .= ' OUTPUT=' . $outputBAM;
    $command .= ' METRICS_FILE=' . $outputMetrics;
    $command .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam($rH_cfg, 'markDup', 'markDupRecInRam', 1, 'int');

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub collectMetrics {
  my $rH_cfg        = shift;
  my $inputBAM      = shift;
  my $outputMetrics = shift;
  my $reference     = shift;

  my $rO_job = new Job();
  $rO_job->testInputOutputs([$inputBAM], [$outputMetrics . '.quality_by_cycle.pdf']);

  if (!$rO_job->isUp2Date()) {
    if (!defined($reference)) {
      $reference = LoadConfig::getParam($rH_cfg, 'collectMetrics', 'referenceFasta', 1, 'filepath');
    }
  
    my $command;
    $rO_job->addModules($rH_cfg, [
      ['collectMetrics', 'moduleVersion.java'],
      ['collectMetrics', 'moduleVersion.picard'],
      ['collectMetrics', 'moduleVersion.cranR']
    ]); 
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'collectMetrics', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'collectMetrics', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'collectMetrics', 'collectMetricsRam') . ' -jar \${PICARD_HOME}/CollectMultipleMetrics.jar';
    $command .= ' PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics  VALIDATION_STRINGENCY=SILENT';
    $command .= ' TMP_DIR=' . LoadConfig::getParam($rH_cfg, 'collectMetrics', 'tmpDir');
    $command .= ' REFERENCE_SEQUENCE=' . $reference;
    $command .= ' INPUT=' . $inputBAM;
    $command .= ' OUTPUT=' . $outputMetrics;
    $command .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam($rH_cfg, 'collectMetrics', 'collectMetricsRecInRam', 1, 'int');

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

# Sort BAM/SAM files
sub sortSam {
  my $rH_cfg    = shift;
  my $input     = shift;
  my $output    = shift;
  my $sortOrder = shift;

  # If command input/output are stdin/stdout, they are ignored
  my $rO_job = new Job($input eq "/dev/stdin" ? [] : [$input], $output eq "/dev/stdout" ? [] : [$output]);
#  $rO_job->testInputOutputs([$input], [$output]);

  if (!$rO_job->isUp2Date()) {
    $rO_job->addModules($rH_cfg, [['sortSam', 'moduleVersion.java'], ['sortSam', 'moduleVersion.picard']]);
    my $command .= "java -Djava.io.tmpdir=" . LoadConfig::getParam($rH_cfg, 'sortSam', 'tmpDir') . " " . LoadConfig::getParam($rH_cfg, 'sortSam', 'extraJavaFlags') . " -Xmx" . LoadConfig::getParam($rH_cfg, 'sortSam', 'sortRam') . " -jar \\\${PICARD_HOME}/SortSam.jar";
    $command .= " \\\n  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true";
    $command .= " \\\n  TMP_DIR=" . LoadConfig::getParam($rH_cfg, 'sortSam', 'tmpDir');
    $command .= " \\\n  INPUT=$input";
    $command .= " \\\n  OUTPUT=$output";
    $command .= " \\\n  SORT_ORDER=$sortOrder";
    $command .= " \\\n  MAX_RECORDS_IN_RAM=" . LoadConfig::getParam($rH_cfg, 'sortSam', 'sortRecInRam', 1, 'int');

    $rO_job->addCommand($command);
  }
  return $rO_job;
}


# reorder BAM/SAM files based on reference/dictionary
sub reorderSam {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBAM      = shift;
  my $outputBAM     = shift;

  my $rO_job = new Job();
  $rO_job->testInputOutputs([$inputBAM], [$outputBAM]);

  if (!$rO_job->isUp2Date()) {
    my $command;
    $rO_job->addModules($rH_cfg, [['reorderSam', 'moduleVersion.java'], ['reorderSam', 'moduleVersion.picard']]);
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'reorderSam', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'reorderSam', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'reorderSam', 'reorderRam') . ' -jar \${PICARD_HOME}/ReorderSam.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true';
    $command .= ' TMP_DIR=' . LoadConfig::getParam($rH_cfg, 'reorderSam', 'tmpDir');
    $command .= ' INPUT=' . $inputBAM;
    $command .= ' OUTPUT=' . $outputBAM;
    $command .= ' REFERENCE=' . LoadConfig::getParam($rH_cfg, 'reorderSam', 'referenceFasta', 1, 'filepath');
    $command .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam($rH_cfg, 'reorderSam', 'reorderRecInRam', 1, 'int');

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

# Convert SAM/BAM file to fastq format
sub samToFastq {
  my $rH_cfg       = shift;
  my $inputSAMBAM  = shift;
  my $outputFastq1 = shift;    # If single end reads, $outputFastq1 will store the FASTQ output
  my $outputFastq2 = shift;    # Used for paired end reads only

  my $rO_job;

  # If paired end reads
  if ($outputFastq2) {
    $rO_job = Job->new([$inputSAMBAM], [$outputFastq1, $outputFastq2]);
  # else single end reads
  } else {
    $rO_job = Job->new([$inputSAMBAM], [$outputFastq1]);
  }

  if (!$rO_job->isUp2Date()) {
    my $command;
    $rO_job->addModules($rH_cfg, [
      ['samToFastq', 'moduleVersion.java'],
      ['samToFastq', 'moduleVersion.picard']
    ]);
    $command .= "java -Djava.io.tmpdir=" . LoadConfig::getParam($rH_cfg, 'samToFastq', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'samToFastq', 'extraJavaFlags') . " -jar \\\${PICARD_HOME}/SamToFastq.jar \\\n";
    $command .= "  INPUT=$inputSAMBAM \\\n";
    $command .= "  FASTQ=$outputFastq1";

    # If paired end reads
    if (defined($outputFastq2)) {
      $command .= " \\\n  SECOND_END_FASTQ=$outputFastq2";
    }

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

1;
