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

package Tools;

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
use Cwd 'abs_path';
use File::Basename;
use LoadConfig;

# SUB
#-----------------------
sub getToolShedDir {
  my $currentDir = dirname(__FILE__);
  return abs_path($currentDir.'/../tool_shed');
}

sub filterNStretches {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $inputVCF    = shift;
  my $outputVCF   = shift;

  my $toolShedDir = getToolShedDir();

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputVCF], [$outputVCF]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['metrics' , 'moduleVersion.tools']]) . ' &&';
    $command .= ' \$PERL_TOOLS/filterLongIndel.pl ';
    $command .= ' ' . $inputVCF;
    $command .= ' > ' . $outputVCF;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

1;
