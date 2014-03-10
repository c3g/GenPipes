#!/usr/bin/env perl

=head1 NAME

I<Celera>

=head1 SYNOPSIS

Celera->run()

=head1 DESCRIPTION

B<Celera> This a library to run Celera assembler. Initially written for PacBio Assembly pipeline.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debug

=cut

package Celera;

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
sub runCelera {
  my $rH_cfg = shift;
  my $outdir = shift;
  my $prefix = shift;
  my $ini    = shift;
  my $infile = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outdir . "/9-terminator/$prefix.ctg.fasta"]); 

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= ' rm -rf ' . $outdir . '/* && ';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [['memtime', 'moduleVersion.memtime'], ['celera', 'moduleVersion.celera']]) . ' && ';
    $cmd .= ' memtime';
    $cmd .= ' runCA';
    $cmd .= ' -d ' . $outdir;
    $cmd .= ' -p ' . $prefix;
    $cmd .= ' -s ' . $ini;
    $cmd .= ' ' . $infile;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub fastqToCA {
  my $rH_cfg      = shift;  
  my $libraryName = shift;
  my $reads       = shift;
  my $outfile     = shift;
  
  my $ro_job = new Job();
  $ro_job->testInputOutputs([$reads], [$outfile]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [['memtime', 'moduleVersion.memtime'], ['celera', 'moduleVersion.celera']]) . ' && ';
    $cmd .= ' memtime';
    $cmd .= ' fastqToCA';
    $cmd .= ' -technology sanger';
    $cmd .= ' -type sanger';
    $cmd .= ' -libraryname ' . $libraryName;
    $cmd .= ' -reads ' . $reads;
    $cmd .= ' > ' . $outfile;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

1;
