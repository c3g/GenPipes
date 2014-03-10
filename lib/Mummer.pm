#!/usr/bin/env perl

=head1 NAME

I<Mummer>

=head1 SYNOPSIS

Mummer->nucmer()

=head1 DESCRIPTION

B<Mummer> This a library to run MUMmer software.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package Mummer;

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
use Job;
use LoadConfig;

# SUB
#-----------------------
sub nucmer {
  my $rH_cfg         = shift;
  my $c              = shift;
  my $prefix         = shift;
  my $fastaRef       = shift;
  my $fastaConsensus = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([""], [""]);

  if (!$ro_job->isUp2Date()) {
    my $cmd;

    $cmd = LoadConfig::moduleLoad($rH_cfg, [['mummer', 'moduleVersion.mummer']]) . ' && ';
    $cmd .= ' nucmer';
    $cmd .= ' -maxmatch';
    $cmd .= ' -c ' . $c;
    $cmd .= ' -p ' . $prefix;
    $cmd .= ' ' . $fastaRef;
    $cmd .= ' ' . $fastaConsensus;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub mummerPlot {
  my $rH_cfg     = shift;
  my $title      = shift;
  my $prefix     = shift;
  my $outfile    = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([""], [""]);

  if (!$ro_job->isUp2Date()) {
    my $cmd;

    $cmd = LoadConfig::moduleLoad($rH_cfg, [['mummer', 'moduleVersion.mummer']]) . ' && ';
    $cmd .= ' mummerplot';
    $cmd .= ' -title ' . $title;
    $cmd .= ' --png';
    $cmd .= ' --layout';  
    $cmd .= ' --filter';
    $cmd .= ' -p ' . $prefix;
    $cmd .= ' ' . $outfile;  

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub dnadiff {
  my $rH_cfg         = shift;
  my $prefix         = shift;
  my $referenceFasta = shift;
  my $consensusFasta = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([""], [""]);

  if (!$ro_job->isUp2Date()) {
    my $cmd;

    $cmd = LoadConfig::moduleLoad($rH_cfg, [['mummer', 'moduleVersion.mummer']]) . ' && ';
    $cmd .= 'dnadiff';
    $cmd .= ' -p ' . $prefix;
    $cmd .= ' ' . $referenceFasta;
    $cmd .= ' ' . $consensusFasta;  

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub showsnp {
  my $rH_cfg  = shift;
  my $x       = shift;
  my $infile  = shift;
  my $outfile = shift;

    my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfile]);

  if (!$ro_job->isUp2Date()) {
    my $cmd;

    $cmd = LoadConfig::moduleLoad($rH_cfg, [['mummer', 'moduleVersion.mummer']]) . ' && ';
    $cmd .= 'show-snps';
    $cmd .= ' -rlTC';
    $cmd .= ' -x ' . $x;
    $cmd .= ' ' . $infile;
    $cmd .= ' > ' . $outfile;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub reference {
  my $rH_cfg         = shift;
  my $prefix1        = shift;
  my $fastaRef       = shift;
  my $fastaConsensus = shift;
  my $title          = shift;
  my $prefix2        = shift;
  my $outfile        = shift;
  my $prefix3        = shift;
  my $infile2        = shift;
  my $outfile2       = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$fastaConsensus, $fastaRef], [$outfile2]);

  if (!$ro_job->isUp2Date()) {
    my $cmd;

    $cmd = LoadConfig::moduleLoad($rH_cfg, [['mummer', 'moduleVersion.mummer']]) . ' &&';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [['gnuplot', 'moduleVersion.gnuplot']]) . ' &&';
    $cmd .= ' promer';
    $cmd .= ' -maxmatch';
    $cmd .= ' -c ' . LoadConfig::getParam($rH_cfg, 'mummer', 'c');
    $cmd .= ' -p ' . $prefix1;
    $cmd .= ' ' . $fastaRef;
    $cmd .= ' ' . $fastaConsensus;
    $cmd .= ' && ';  
    $cmd .= ' mummerplot';
    $cmd .= ' -title ' . $title;
    $cmd .= ' --png';
    $cmd .= ' --layout';  
    $cmd .= ' --filter';
    $cmd .= ' -p ' . $prefix2;
    $cmd .= ' ' . $outfile;  
    $cmd .= ' && ';
    $cmd .= 'dnadiff';
    $cmd .= ' -p ' . $prefix3;
    $cmd .= ' ' . $fastaRef;
    $cmd .= ' ' . $fastaConsensus;  
    $cmd .= ' && ';
    $cmd .= 'show-snps';
    $cmd .= ' -rlTC';
    $cmd .= ' -x ' . LoadConfig::getParam($rH_cfg, 'mummer', 'x');
    $cmd .= ' ' . $infile2;
    $cmd .= ' > ' . $outfile2;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub self {
  my $rH_cfg         = shift;
  my $prefix1        = shift;
  my $fastaConsensus = shift;
  my $title          = shift;
  my $prefix2        = shift;
  my $outfile        = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$fastaConsensus], [$outfile]);

  if (!$ro_job->isUp2Date()) {
    my $cmd;

    $cmd = LoadConfig::moduleLoad($rH_cfg, [['mummer', 'moduleVersion.mummer']]) . ' &&';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [['gnuplot', 'moduleVersion.gnuplot']]) . ' &&';
    $cmd .= ' nucmer';
    $cmd .= ' -maxmatch';
    $cmd .= ' -c ' . LoadConfig::getParam($rH_cfg, 'mummer', 'c');
    $cmd .= ' -p ' . $prefix1;
    $cmd .= ' ' . $fastaConsensus;
    $cmd .= ' ' . $fastaConsensus;
    $cmd .= ' && ';  
    $cmd .= ' mummerplot';
    $cmd .= ' -title ' . $title;
    $cmd .= ' --png';
    $cmd .= ' --layout';  
    $cmd .= ' --filter';
    $cmd .= ' -p ' . $prefix2;
    $cmd .= ' ' . $outfile;  

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

1;
