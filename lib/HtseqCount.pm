#!/usr/bin/env perl

=head1 NAME

I<HtseqCount>

=head1 SYNOPSIS

HtseqCount:sub(args)

B<HtseqCount::readCount>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group) 

B<HtseqCount::sortRead>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group)

B<HtseqCount::matrixMake>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group)

All subroutines return a ref_hash with the command line



=head1 DESCRIPTION

B<HtseqCount> is a library to generate basic statistics on raw read count.

=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk
Mathieu Bourgey mbourgey@genomequebec.com

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug


=cut

package HtseqCount;

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
use Config::Simple;
use Data::Dumper;
use File::Basename;
use LoadConfig;
use SAMtools;

our $rH_cfg;
our $sampleName;
our $rH_laneInfo;
our $readFile;
our $group;

sub matrixMake {
  $rH_cfg      = shift;
  $sampleName  = shift;
  $rH_laneInfo = shift;
  my $db       = shift;    # blast database
  $group       = shift;

  $group = ( !defined $group ) ? $sampleName : $group;


  # option used if more than one db was specified on the config file.
  # In this case $db should be passed as an argument
  #------------------------------------------------------------------
  $rH_cfg->{'blast.db'} = defined($db) ? $db : $rH_cfg->{'blast.db'};

  my %retVal;
  my $command       = '';
  my $laneDirectory = "read_count/" . $group . "/";
  my $input = 'assembly/' . $group . '/fasta_split/' . basename($rH_cfg->{'blast.db'}) . '/blast_BestHit.txt';

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$input], [$laneDirectory . 'tmpmatrix.csv DGE/' . $group . '/matrix.csv']);

  if (!$ro_job->isUp2Date()) {
    $command .= ' sh ' . $rH_cfg->{'htseq.tempMatrix'} . ' ' . 'alignment/' . $group . '/' . $group . '.gtf';
    $command .= ' ' . $input;
    $command .= ' ' . $laneDirectory . 'tmpmatrix.csv &&';
    $command .= ' sh ' . $rH_cfg->{'htseq.fullMatrix'} . ' ' . $laneDirectory . ' &&';
    $command .= ' cp ' . $laneDirectory . 'tmpmatrix.csv DGE/' . $group . '/matrix.csv ';

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub readCountPortable {
  my $rH_cfg     = shift;
  my $inputBam   = shift;
  my $inputGtf   = shift;
  my $outputFile = shift;
  my $strandInfo = shift;

  if (!(defined $strandInfo)) {
    $strandInfo='no';
  }

  my $command;
  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputBam], [$outputFile]);

  if (!$ro_job->isUp2Date()) {
    #Just need the command
    my $rO_viewFilterJob = SAMtools::viewFilter($rH_cfg, $inputBam);

    $command .= LoadConfig::moduleLoad($rH_cfg, [['htseq', 'moduleVersion.python']]) . ' && ';
    $command .= ' ' . $rO_viewFilterJob->getCommand(0);
    $command .= ' | htseq-count - ' .  $inputGtf;
    $command .= ' -s ' . $strandInfo;
    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'htseq','options');
    $command .= ' >' . $outputFile;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}


sub refGtf2matrix {
  my $rH_cfg      = shift;
  my $refGtf  = shift;
  my $readCountDir = shift;
  my $readcountExtension = shift;
  my $outputDir = shift;
  my $outputMatrix  = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$refGtf], [$outputDir . '/' . $outputMatrix]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['htseq', 'moduleVersion.tools']]) . ' &&';
    $command .= ' gtf2tmpMatrix.awk ' . $refGtf;
    $command .= ' ' . $outputDir . '/tmpMatrix.txt &&';
    $command .= ' HEAD=\"Gene\tSymbol\" &&';
    $command .= ' for i in \` ls ' . $readCountDir . '/*' . $readcountExtension . ' \` ;';
    $command .= ' do sort -k1,1 \$i > ' . $outputDir . '/tmpSort.txt ;';
    $command .= ' join -1 1 -2 1 ' . $outputDir . '/tmpMatrix.txt ' . $outputDir . '/tmpSort.txt > ' . $outputDir . '/tmpMatrix.2.txt ;';
    $command .= ' mv ' . $outputDir . '/tmpMatrix.2.txt ' . $outputDir . '/tmpMatrix.txt ;';
    $command .= ' na=\$(basename \$i ' . $readcountExtension . ') ;'; 
    $command .= ' HEAD=\"\$HEAD\t\$na\" ;';
    $command .= ' done &&';
    $command .= ' echo -e \$HEAD | cat - ' . $outputDir . '/tmpMatrix.txt | tr \' \' \'\t\' > ' . $outputDir . '/' . $outputMatrix . ' &&';
    $command .= ' rm ' . $outputDir . '/tmpSort.txt ' . $outputDir . '/tmpMatrix.txt ';

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

1;
