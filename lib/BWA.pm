#!/usr/bin/env perl

=head1 NAME

I<BWA>

=head1 SYNOPSIS

BWA->aln()

=head1 DESCRIPTION

B<BWA> is a library that aligns fastqs on a reference genome

Input = file_name

Output = array


=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package BWA;

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
use LoadConfig;

# SUB
#-----------------------

sub mem {
  my $rH_cfg    = shift;
  my $idxBase = shift;
  my $in1Fastq  = shift;
  my $in2Fastq  = shift;
  my $outSam = shift;
  my $readGroup = shift;

  my $rO_job = new Job([$idxBase, $in1Fastq, $in2Fastq], [$outSam]);

  my $command = "";
  if ($outSam) {
    # Create output directory
    $command = "mkdir -p " . dirname($outSam) . " && \\\n";
  }

  $rO_job->addModules($rH_cfg, [['mem', 'moduleVersion.bwa']]);
  $command .= "bwa mem";
  $command .= " \\\n  -t " . LoadConfig::getParam($rH_cfg, 'mem', 'bwaMemThreads', 1, 'int');
  $command .= " \\\n  " . LoadConfig::getParam($rH_cfg, 'mem', 'bwaExtraFlags');
  if ($readGroup) {
    $command .= " \\\n  -R " . $readGroup;
  }
  $command .= " \\\n  " . $idxBase;
  $command .= " \\\n  " . $in1Fastq;
  if ($in2Fastq) {
    $command .= " \\\n  " . $in2Fastq;
  }
  if ($outSam) {
    $command .= " \\\n  > " . $outSam;
  }
  $rO_job->addCommand($command);

  return $rO_job;
}

sub aln {
  my $rH_cfg       = shift;
  my $inDbFasta    = shift;
  my $inQueryFastq = shift;
  my $outSai       = shift;

  my $rO_job = new Job([$inDbFasta, $inQueryFastq], [$outSai]);

  # Create output directory
  my $command = "mkdir -p " . dirname($outSai) . " && \\\n";

  $rO_job->addModules($rH_cfg, [['aln', 'moduleVersion.bwa']]);
  $command .= "bwa aln";
  $command .= " \\\n  -t " . LoadConfig::getParam($rH_cfg, 'aln', 'bwaAlnThreads', 1, 'int');
  $command .= " \\\n  " . $inDbFasta;
  $command .= " \\\n  " . $inQueryFastq;
  if ($outSai) {
    $command .= " \\\n  -f " . $outSai;
  }
  $rO_job->addCommand($command);

  return $rO_job;
}

sub sampe {
  my $rH_cfg    = shift;
  my $inDbFasta = shift;
  my $in1Sai    = shift;
  my $in2Sai    = shift;
  my $in1Fastq  = shift;
  my $in2Fastq  = shift;
  my $outSam    = shift;
  my $readGroup = shift;

  my $rO_job = new Job([$inDbFasta, $in1Sai, $in2Sai, $in1Fastq, $in2Fastq], [$outSam]);

  $rO_job->addModules($rH_cfg, [['aln', 'moduleVersion.bwa']]);
  my $command .= "bwa sampe";
  $command .= " " . LoadConfig::getParam($rH_cfg, 'aln', 'bwaExtraSamXeFlags', 0);
  if ($readGroup) {
    $command .= " \\\n  -r " . $readGroup;
  }
  $command .= " \\\n  " . $inDbFasta;
  $command .= " \\\n  " . $in1Sai;
  $command .= " \\\n  " . $in2Sai;
  $command .= " \\\n  " . $in1Fastq;
  $command .= " \\\n  " . $in2Fastq;
  if ($outSam) {
    $command .= " \\\n  -f " . $outSam;
  }

  $rO_job->addCommand($command);

  return $rO_job;
}

sub samse {
  my $rH_cfg    = shift;
  my $inDbFasta = shift;
  my $inSai     = shift;
  my $inFastq   = shift;
  my $outSam    = shift;
  my $readGroup = shift;

  my $rO_job = new Job([$inDbFasta, $inSai, $inFastq], [$outSam]);

  $rO_job->addModules($rH_cfg, [['aln', 'moduleVersion.bwa']]);
  my $command .= "bwa samse ";
  $command .= " " . LoadConfig::getParam($rH_cfg, 'aln', 'bwaExtraSamXeFlags', 0);
  if ($readGroup) {
    $command .= " \\\n  -r " . $readGroup;
  }
  $command .= " \\\n  " . $inDbFasta;
  $command .= " \\\n  " . $inSai;
  $command .= " \\\n  " . $inFastq;
  if ($outSam) {
    $command .= " \\\n  -f " . $outSam;
  }

  $rO_job->addCommand($command);

  return $rO_job;
}

1;
