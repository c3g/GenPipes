#!/usr/env/perl

=head1 NAME

I<BLAST>

=head1 SYNOPSIS

B<BLAST::align>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $FastaFile)


B<BLAST::alignParallel>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $FastaFile)

B<BLAST::alignParallel> will depend on two scripts and their paths must be put on the config file

1- B<ParallelBlast>: a script written by David Morais 

2- B<fastasplit>: www.bcgsc.ca/downloads/parts/software/resources/src/exonerate-1.4.0/src/util/fastasplit.c

B<Both scripts need to be in the same dierctory.>

=head1 DESCRIPTION

B<BLAST> is a library to use the
alignment package, BLAST.

The lib implements two subroutines: one that performs BLAST with its on parallelization options;
and another that performs an external parallelization. 

In my benchmark the external parallelization was 25% faster but fell free to try it out yourself.


=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug




=cut

package BLAST;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Data::Dumper;
use Config::Simple;
use LoadConfig;

#-------------------
# SUB
#-------------------
our $rH_cfg;
our $sampleName;
our $rH_laneInfo;
our $fileFasta;

sub align {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $fileFasta   = shift;

    my $outFile = $fileFasta;
    $outFile =~ s/\*_//;
    my $command = '';
    my %retVal;
    my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";

    $command .= 'module add mugqic/BLAST/2.2.26+;';
    $command .= ' ' . $rH_cfg->{'blast.program'} . ' -num_threads ' . $rH_cfg->{'blast.nbThreads'};
    $command .= ' -query ' . $laneDirectory . 'fasta_split/' . $fileFasta;
    $command .= ' -db ' . $rH_cfg->{'blast.db'} . ' -out ' . $laneDirectory . 'fasta_split/' . $outFile . '_BLASTOUT.txt';
    $command .= ' ' . $rH_cfg->{'blast.options'};

    $retVal{'command'} = $command;
    return ( \%retVal );

}

sub alignParallel {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $fileFasta   = shift;

    my $outFile = $fileFasta;
    $outFile =~ s/\*_//;
    my $command = '';
    my %retVal;
    my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";

    $command .= 'module add mugqic/BLAST/2.2.26+;';
    $command .= ' ' . $rH_cfg->{'blast.parallelBlast'} . ' --file ' . $laneDirectory . 'fasta_split/' . $fileFasta;
    $command .= ' --OUT ' . $laneDirectory . 'fasta_split/' . $outFile . '_BLASTOUT.txt';
    $command .= ' -n ' . $rH_cfg->{'blast.nbThreads'} . ' --BLAST ' . '\'' . $rH_cfg->{'blast.program'};
    $command .= ' -db ' . $rH_cfg->{'blast.db'} . ' ' . $rH_cfg->{'blast.options'} . '\'';

    $retVal{'command'} = $command;
    return ( \%retVal );

}

1;

