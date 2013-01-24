#!/usr/env/perl

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

# Dependencies
#-----------------------
use Data::Dumper;
use Config::Simple;
use File::Basename;
our $rH_cfg;
our $sampleName;
our $rH_laneInfo;
our $readFile;
our $group;

sub matrixMake {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    my $db = shift;    # blast database
    $group       = shift;
    
	$group = ( !defined $group ) ? $sampleName : $group;
    

    # option used if more than one db was specified on the config file.
    # In this case $db should be passed as an argument
    #------------------------------------------------------------------
    $rH_cfg->{'blast.db'} = defined($db) ? $db : $rH_cfg->{'blast.db'};

    my %retVal;
    my $command       = '';
    my $laneDirectory = "read_count/" . $group . "/";

    $command .= ' sh ' . $rH_cfg->{'htseq.tempMatrix'} . ' ' . 'alignment/' . $group . '/' . $group . '.gtf';
    $command .= ' ' . 'assembly/' . $group . '/fasta_split/' . basename($rH_cfg->{'blast.db'}) . '/blast_BestHit.txt ';
    $command .= ' ' . $laneDirectory . 'tmpmatrix.csv ;';
    $command .= ' sh ' . $rH_cfg->{'htseq.fullMatrix'} . ' ' . $laneDirectory . ' ;';
    $command .= ' cp ' . $laneDirectory . 'tmpmatrix.csv DGE/' . $group . '/matrix.csv ;';

    $retVal{'command'} = $command;
    return ( \%retVal );
}

sub readCount {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $group       = shift;
	$group = ( !defined $group ) ? $sampleName : $group;
	
    my %retVal;
    my $command       = '';
    my $laneDirectory = "read_count/" . $group . "/";

    $command .= ' module add mugqic/samtools/0.1.6; ';
    $command .= ' samtools view ' . $laneDirectory . $sampleName . '.QueryName.bam | ';
    $command .= ' htseq-count - ' . 'alignment/' . $group . '/' . $group . '.gtf ';
    $command .= ' -s no >' . $laneDirectory . $sampleName . '.readcount.cvs';

    $retVal{'command'} = $command;
    return ( \%retVal );
}

sub sortRead {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $group       = shift;
	
	$group = ( !defined $group ) ? $sampleName : $group;
    my %retVal;
    my $command       = '';
    my $laneDirectory = "alignment/" . $group . "/";

    $command .= ' mkdir -p  read_count/' . $group . ' ;';
    $command .= ' mkdir -p  DGE/' . $group . ';';
    $command .= ' module add jdk64 ; module add mugqic/picard/1.64/ ; ';
    $command .= ' java -Xmx30g -jar SortSam.jar ';
    $command .= ' VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true TMP_DIR=/mnt/scratch_mp2/bourque/bourque_group/tmpDir/ ';
    $command .= ' INPUT=' . $laneDirectory . $sampleName . '.sorted.bam ';
    $command .= ' OUTPUT=read_count/' . $group . '/' . $sampleName . '.QueryName.bam ';
    $command .= ' SORT_ORDER=queryname ';

    $retVal{'command'} = $command;
    return ( \%retVal );
}
1;

