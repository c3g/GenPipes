#!/usr/env/perl

=head1 NAME

I<Pindel>

=head1 SYNOPSIS

Pindel

=head1 DESCRIPTION

B<Cfreec> is a library to analyse CNV events in a genome

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package Cfreec;

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
sub pairedFreec {
  my $rH_cfg       = shift;
  my $tumorPileup  = shift;
  my $normalPileup = shift;
  my $sampleConfig = shift;
  my $output       = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$normalPileup, $tumorPileup], [$sampleConfig]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['ControlFreec', 'moduleVersion.controlFreec']]) . ' && ';
    $command .= 'sed \"s|TUMOR_PILEUP|' . $tumorPileup . '|g\" ' . '\${FREEC_HOME}/' . LoadConfig::getParam($rH_cfg, 'ControlFreec', 'referenceConfigFile') . ' > ' . $sampleConfig . ' && ';
    $command .= 'sed \"s|NORMAL_PILEUP|' . $normalPileup . '|g\" -i ' . $sampleConfig . ' && ';
    $command .= 'sed \"s|OUTPUT_DIR|' . $output . '|g\"  -i ' . $sampleConfig . ' && ';
    $command .= 'sed \"s|FORMAT_TYPE|' . LoadConfig::getParam($rH_cfg, 'ControlFreec', 'inputType') . '|g\"  -i ' . $sampleConfig . ' && ';
    $command .= 'freec -conf ' . $sampleConfig ;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

1;
