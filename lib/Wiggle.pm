#!/usr/bin/env perl

=head1 NAME

I<Wiggle>

=head1 SYNOPSIS

Wiggle-> strandBam()
Wiggle-> graph()

=head1 DESCRIPTION

B<Wiggle> is a library to generate wiggle graphs

Input = file_name

Output = array


=head1 AUTHOR
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package Wiggle;

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
use Picard;

# SUB
#-----------------------
sub strandBam {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBAM      = shift;
  my $rA_outputBAM  = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputBAM], [$rA_outputBAM->[0], $rA_outputBAM->[1]]);

  my @mergeBAMFtmp = ($inputBAM . 'tmp1.forward.bam', $inputBAM . 'tmp2.forward.bam');
  my @mergeBAMRtmp = ($inputBAM . 'tmp1.reverse.bam', $inputBAM . 'tmp2.reverse.bam');
  my $mergeFJob = Picard::mergeFiles($rH_cfg, $sampleName, \@mergeBAMFtmp, $rA_outputBAM->[0]);
  my $mergeRJob = Picard::mergeFiles($rH_cfg, $sampleName, \@mergeBAMRtmp, $rA_outputBAM->[1]);

  if (!$ro_job->isUp2Date()) {
    my $Fcommand = LoadConfig::moduleLoad($rH_cfg, [['wiggle', 'moduleVersion.samtools']]) . ' &&';
    $Fcommand .= ' samtools view -bh -F 256 -f 81 ' . $inputBAM;
    $Fcommand .= ' > ' . $inputBAM . 'tmp1.forward.bam &&';
    $Fcommand .= ' samtools view -bh -F 256 -f 161 ' . $inputBAM;
    $Fcommand .= ' > ' . $inputBAM . 'tmp2.forward.bam &&';
    $Fcommand .= ' ' . $mergeFJob->getCommand(0) . ' && ';
    $Fcommand .= ' rm ' . $inputBAM . 'tmp*.forward.*am';

    $ro_job->addCommand($Fcommand);

    my $Rcommand = LoadConfig::moduleLoad($rH_cfg, [['wiggle', 'moduleVersion.samtools']]) . ' &&';
    $Rcommand .= ' samtools view -bh -F 256 -f 97 ' . $inputBAM;
    $Rcommand .= ' > ' . $inputBAM . 'tmp1.reverse.bam &&';
    $Rcommand .= ' samtools view -bh -F 256 -f 145 ' . $inputBAM;
    $Rcommand .= ' > ' . $inputBAM . 'tmp2.reverse.bam &&';
    $Rcommand .= ' ' . $mergeRJob->getCommand(0) . ' &&';
    $Rcommand .= ' rm ' . $inputBAM . 'tmp*.reverse.*am';

    $ro_job->addCommand($Rcommand);
  }

  return $ro_job;
}

sub graph {
  my $rH_cfg         = shift;
  my $sampleName     = shift;
  my $inputBAM       = shift;
  my $outputBegGraph = shift;
  my $outputWiggle   = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputBAM], [$outputBegGraph, $outputWiggle]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [
      ['wiggle', 'moduleVersion.samtools'],
      ['wiggle', 'moduleVersion.bedtools'],
      ['wiggle', 'moduleVersion.bed2wig']
    ]) . ' &&';
    $command .= ' nmblines=\$(samtools view -F 256 -f 81 ' . $inputBAM . ' | wc -l) &&';
    $command .= ' scalefactor=0\$(echo \"scale=2; 1 / (\$nmblines / 10000000);\" | bc) &&';
    $command .= ' genomeCoverageBed -bg -ibam ' . $inputBAM;
    $command .= ' -g ' . LoadConfig::getParam($rH_cfg, 'wiggle', 'chromosomeSizeFile', 1, 'filepath');
    $command .= ' -split -scale \$scalefactor ';
    $command .= ' > ' . $outputBegGraph . ' &&';
    $command .= ' bedGraphToBigWig ' . $outputBegGraph;
    $command .= '  ' . LoadConfig::getParam($rH_cfg, 'wiggle', 'chromosomeSizeFile', 1, 'filepath');
    $command .= '  ' . $outputWiggle;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

sub zipWig {
  my $rH_cfg        = shift;
  my $wigFolder     = shift;
  my $wigArchive    = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs(undef, undef);

  if (!$ro_job->isUp2Date()) {
    my $command = ' zip -r ' . $wigArchive . ' ' . $wigFolder;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

1;
