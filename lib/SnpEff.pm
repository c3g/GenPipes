#!/usr/env/perl

=head1 NAME

I<SVTools>

=head1 SYNOPSIS

Picard->merge()

=head1 DESCRIPTION

B<SVTools> is a library to analyse BAMs for Structural Variants

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package SnpEff;

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
sub annotateDbSnp {
  my $rH_cfg      = shift;
  my $inputVCF    = shift;
  my $outputVCF   = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputVCF], [$outputVCF]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['annotateDbSnp', 'moduleVersion.java'], ['annotateDbSnp', 'moduleVersion.snpeff']]) . ' &&';
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'siftRam');
    $command .= ' -jar \${SNPEFF_HOME}/SnpSift.jar annotate';
    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'dbSnp');
    $command .= ' ' . $inputVCF;
    $command .= ' > ' . $outputVCF;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub computeEffects {
  my $rH_cfg            = shift;
  my $inputVCF          = shift;
  my $outputVCF         = shift;
  my $addEffectSplitter = shift;

  my $addSplit = 0;
  if (defined($addEffectSplitter) && $addEffectSplitter == 1) {
    $addSplit = 1;
  }

  my $ro_job = new Job();
  if ($addSplit == 1) {
    $ro_job->testInputOutputs([$inputVCF], [$outputVCF, $outputVCF . '.statsFile.txt']);
  } else {
    $ro_job->testInputOutputs([$inputVCF], [$outputVCF]);
  }

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['computeEffects', 'moduleVersion.java'], ['computeEffects', 'moduleVersion.snpeff']]);
    if ($addSplit == 1) {
      $command .= ' ' . LoadConfig::getParam($rH_cfg, 'computeEffects', 'moduleVersion.tools');
    }
    $command .= ' &&';
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'computeEffects', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'computeEffects', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'computeEffects', 'snpeffRam');
    $command .= ' -jar \${SNPEFF_HOME}/snpEff.jar eff';
    $command .= ' -c \${SNPEFF_HOME}/snpEff.config';
    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'computeEffects', 'snpeffParams');
    $command .= ' -o vcf';
    $command .= ' -i vcf';
    $command .= ' -csvStats';
    $command .= ' -stats ' . $outputVCF . '.stats.csv';
    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'computeEffects', 'referenceSnpEffGenome');
    $command .= ' ' . $inputVCF;
    $command .= ' > ' . $outputVCF;
    if ($addSplit == 1) {
      $command .= ' && splitSnpEffStat.awk ' . $outputVCF . '.stats.csv';
      $command .= ' ' . $outputVCF . '.part ' . $outputVCF . '.statsFile.txt';
    }

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub annotateDbNSFP {
  my $rH_cfg      = shift;
  my $inputVCF    = shift;
  my $outputVCF   = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputVCF], [$outputVCF]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['annotateDbNSFP', 'moduleVersion.java'], ['annotateDbNSFP', 'moduleVersion.snpeff']]) . ' &&';
    $command .= ' java -Djava.io.tmpdir=' . LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'tmpDir') . ' ' . LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'extraJavaFlags') . ' -Xmx' . LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'siftRam');
    $command .= ' -jar \${SNPEFF_HOME}/SnpSift.jar dbnsfp';
    $command .= ' -v ' . LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'dbNSFP');
    $command .= ' ' . $inputVCF;
    $command .= ' > ' . $outputVCF;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

1;
