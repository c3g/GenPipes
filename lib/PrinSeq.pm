#!/usr/bin/env perl

=head1 NAME

I<PrinSeq>

=head1 SYNOPSIS

PrinSeq->fatqToFasta()

=head1 DESCRIPTION

B<PrinSeq> This a library to analyze PAcBio data using the SmrtAnalysis suite.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package PrinSeq;

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
sub fastqToFasta {
  my $rH_cfg     = shift;
  my $infile     = shift;
  my $outfile    = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfile . ".fasta"]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    # Fastq to Fasta
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [['prinseq', 'moduleVersion.memtime'], ['prinseq', 'moduleVersion.prinseq']]) . ' &&';
    $cmd .= ' memtime';
    $cmd .= ' prinseq-lite.pl';
    $cmd .= ' -verbose';
    $cmd .= ' -fastq ' . $infile;
    $cmd .= ' -out_format 1';
    $cmd .= ' -out_good ' . $outfile;
    $cmd .= ' > /dev/null 2>&1';

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

1;
