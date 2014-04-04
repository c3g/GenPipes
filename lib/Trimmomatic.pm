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

sub trim {
  my $rH_cfg = shift;
  my $input1 = shift;   # FASTQ pair1 or FASTQ single
  my $input2 = shift;   # FASTQ pair2 for PAIRED_END reads only, undefined for SINGLE_END reads
  my $pairedOutput1 = shift;
  my $unpairedOutput1 = shift;
  my $pairedOutput2 = shift;
  my $unpairedOutput2 = shift;
  my $singleOutput = shift;
  my $qualityOffset = shift;
  my $trimLog = shift;
  my $trimStats = shift;

  my $inputs;
  my $outputs;

  if ($input2) {    # Paired end reads
    $inputs = [$input1, $input2];
    $outputs = [$pairedOutput1, $unpairedOutput1, $pairedOutput2, $unpairedOutput2];
  } else {    # Single end reads
    $inputs = [$input1];
    $outputs = [$singleOutput];
  }

  push(@$outputs, $trimLog, $trimStats);

#  $rO_job->testInputOutputs($inputs, $outputs);
  my $rO_job = new Job($inputs, $outputs);

  if (!$rO_job->isUp2Date()) {
    # Create output directories (remove duplicates from output directory list if necessary)
    my $command = "mkdir -p " . join(" ", keys(%{{map {$_ => 1} map(dirname($_), @$outputs)}})) . " && \\\n";

    $rO_job->addModules($rH_cfg, [
      ['trim', 'moduleVersion.java'],
      ['trim', 'moduleVersion.trimmomatic']
    ]);

    $command .= "java -XX:ParallelGCThreads=1 -Xmx2G -jar \\\$TRIMMOMATIC_JAR " . ($input2 ? "PE" : "SE") . " \\\n";
    $command .= "  -threads " . LoadConfig::getParam($rH_cfg, 'trim', 'nbThreads', 1, 'int') . " \\\n";
    $command .= "  -phred" . ($qualityOffset eq "64" ? $qualityOffset : "33") . " \\\n";
    $command .= "  " . join(" \\\n  ", @$inputs) . " \\\n";
    $command .= "  " . join(" \\\n  ", ($pairedOutput1, $unpairedOutput1, $pairedOutput2, $unpairedOutput2)) . " \\\n";
    $command .= "  ILLUMINACLIP:" . LoadConfig::getParam($rH_cfg, 'trim', 'adapterFile', 1, 'filepath') . LoadConfig::getParam($rH_cfg, 'trim', 'clipSettings') . " \\\n";
    my $minQuality  = LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int');
    if ($minQuality > 0) {
      $command .= "  TRAILING:$minQuality \\\n";
    }
    my $headcrop    = LoadConfig::getParam($rH_cfg, 'trim', 'headcrop', 0, 'int');
    if ($headcrop and $headcrop > 0) {
      $command .= "  HEADCROP:$headcrop \\\n";
    }
    my $minLength   = LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int');
    $command .= "  MINLEN:$minLength \\\n";
    $command .= $qualityOffset eq "64" ? "  TOPHRED33 \\\n" : "";
    $command .= "  2> $trimLog && \\\n";

    # Compute trim stats
    if ($input2) {    # Paired end reads
      $command .= 'grep \"^Input Read\" ' . $trimLog . '| sed \'s/Input Read Pairs: \\([0-9]\\+\\).*Both Surviving: \\([0-9]\\+\\).*Forward Only Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\3/g\' | tr \'#\' \'\n\' > ' . $trimStats;
    } else {    # Single end reads
      $command .= 'grep \"^Input Read\" ' . $trimLog . '| sed \'s/Input Reads: \\([0-9]\\+\\).*Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\2/g\' | tr \'#\' \'\n\' > ' . $trimStats;
    }

    $rO_job->addCommand($command);
  }

  return $rO_job;
}

1;
