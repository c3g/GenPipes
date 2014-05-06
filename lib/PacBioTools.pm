#!/usr/bin/env perl

=head1 NAME

I<PacBioTools>

=head1 SYNOPSIS

PacBioTools->run()

=head1 DESCRIPTION

B<PacBioTools> This a library to analyze PacBio data using the SmrtAnalysis suite.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debug

=cut

package PacBioTools;

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
sub getCutoff {
  my $rH_cfg           = shift;
  my $infile           = shift;
  my $coverage         = shift;
  my $genomeSize       = shift;
  my $coverageFraction = shift;
  #my $xml              = shift;
  #my $xmlOut           = shift;
  my $outfile          = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfile]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';

    # Choose a subread length threshold such that subreads above the threshold provide about 20x coverage of the genome.
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['default', 'moduleVersion.mugqictools'],
      ['default', 'moduleVersion.perl']
    ]) . ' &&';
    $cmd .= ' memtime';
    $cmd .= ' pacBioGetCutoff.pl';
    $cmd .= ' --infile ' . $infile;
    $cmd .= ' --coverage ' . $coverage;
    $cmd .= ' --genomeSize ' . $genomeSize;
    $cmd .= ' --coverageCutoff';
    $cmd .= ' --coverageFraction ' . $coverageFraction;
    #$cmd .= ' --coverageFraction ' . LoadConfig::getParam($rH_cfg, 'preassembly', 'coverageFraction');
    #$cmd .= ' --xml ' . $xml;
    #$cmd .= ' --xmlOut ' . $xmlOut;
    $cmd .= ' > ' . $outfile;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub celeraConfig {
  my $rH_cfg         = shift;
  my $merSize        = shift;
  my $infile         = shift;
  my $outfile        = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfile]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['default', 'moduleVersion.mugqictools'],
      ['default', 'moduleVersion.perl']
    ]) . ' &&';
    $cmd .= ' memtime';
    $cmd .= ' pacBioAssemblyCeleraConfig.pl';
    $cmd .= ' --infile '             . $infile;
    $cmd .= ' --merylThreads '       . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'merylThreads', 1, 'int');
    $cmd .= ' --ovlThreads '         . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'ovlThreads', 1, 'int');
    $cmd .= ' --overlapper '         . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'overlapper');
    $cmd .= ' --merCompression '     . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'merCompression');
    $cmd .= ' --merSize '            . $merSize;
    $cmd .= ' --merylMemory '        . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'merylMemory');
    $cmd .= ' --ovlErrorRate '       . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'ovlErrorRate', 1, 'float');
    $cmd .= ' --ovlMinLen '          . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'ovlMinLen', 1, 'int');
    $cmd .= ' --frgMinLen '          . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'frgMinLen', 1, 'int');
    $cmd .= ' --ovlStoreMemory '     . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'ovlStoreMemory', 1, 'int');
    $cmd .= ' --ovlConcurrency '     . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'ovlConcurrency');
    $cmd .= ' --ovlCorrConcurrency ' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'ovlCorrConcurrency', 1, 'int');
    $cmd .= ' --cnsConcurrency '     . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'cnsConcurrency', 1, 'int');
    $cmd .= ' --frgCorrThreads '     . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'frgCorrThreads', 1, 'int');
    $cmd .= ' --stopAfter '          . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'stopAfter');
    $cmd .= ' --unitigger '          . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'unitigger');
    $cmd .= ' > ' . $outfile;
    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub assemblyStats {
  my $rH_cfg                = shift;
  my $shortReads            = shift;
  my $longReads             = shift;
  my $correctedReads        = shift;
  my $filteredSummary       = shift;
  #my $assemblyQc            = shift;
  my $contigs               = shift;
  my $sampleName            = shift;
  my $suffix                = shift;
  my $estimatedGenomeSize   = shift;
  my $smrtCells             = shift;
  my $outdir                = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs(
    [$filteredSummary],
    [$outdir . "/summaryTableAssembly.tsv", $outdir . "/summaryTableReads.tsv"]
  );

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['default', 'moduleVersion.R'],
      ['default', 'moduleVersion.mugqictools'],
      ['default', 'moduleVersion.perl']
    ]) . ' &&';
    $cmd .= ' memtime ';
    $cmd .= ' pacBioAssemblyStats.pl';
    $cmd .= ' --shortReads ' . $shortReads;
    $cmd .= ' --longReads ' . $longReads;
    $cmd .= ' --correctedReads ' . $correctedReads;
    $cmd .= ' --filteredSummary ' . $filteredSummary;
    #$cmd .= ' --assemblyQc ' . $assemblyQc;
    $cmd .= ' --contigs ' . $contigs;
    $cmd .= ' --sampleName ' . $sampleName;
    $cmd .= ' --suffix ' . $suffix;
    $cmd .= ' --estimatedGenomeSize ' . $estimatedGenomeSize;
    $cmd .= ' --smrtCells ' . $smrtCells;
    $cmd .= ' --outdir ' . $outdir;

    $ro_job->addCommand($cmd);

  }
  return $ro_job;
}

sub splitReads{
  my $rH_cfg              = shift;
  my $subreads            = shift;
  my $cutoff              = shift; # a file containing the cutoff number is actually passed in arg here.
  my $shortReads          = shift;
  my $longReads           = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$subreads], [$shortReads, $longReads]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['default', 'moduleVersion.mugqictools'],
      ['default', 'moduleVersion.perl']
    ]) . ' &&';
    $cmd .= ' memtime';
    $cmd .= ' pacBioSplitReads.pl';
    $cmd .= ' --infile ' . $subreads;
    $cmd .= ' --cutoff \`cat ' . $cutoff . '\` ';
    $cmd .= ' --outfileShort ' . $shortReads;
    $cmd .= ' --outfileLong ' . $longReads;

    $ro_job->addCommand($cmd);

  }
  return $ro_job;
}

sub compile{
  my $rH_cfg              = shift;
  my $indir               = shift;
  my $sampleName          = shift;
  my $estimatedGenomeSize = shift;
  my $outfile             = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs(undef, undef);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['default', 'moduleVersion.mugqictools'],
      ['default', 'moduleVersion.perl']
    ]) . ' ;';
    $cmd .= ' memtime';
    $cmd .= ' pacBioCompileStats.pl';
    $cmd .= ' --indir ' . $indir;
    $cmd .= ' --estimatedGenomeSize ' . $estimatedGenomeSize;
    $cmd .= ' --sampleName ' . $sampleName;
    $cmd .= ' > ' . $outfile;

    $ro_job->addCommand($cmd);

  }
  return $ro_job;
}

1;
