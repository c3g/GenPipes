#!/usr/env/perl

=head1 NAME

I<RemoveDuplicateReads>

=head1 SYNOPSIS

RemoveDuplicateReads::DeleteDuplicates(%ref_hash_config, $sample_name, %ref_hash_laneInfo)

=head1 DESCRIPTION

B<RemoveDuplicateReads> is a library used by delete
duplicate reads

It returns a reference hashe (%retVal);


=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug



=cut

package RemoveDuplicateReads;

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

sub deleteDuplicates {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;

    my $rH_retVal;
    if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
        $rH_retVal = _singleCommand();
    }
    elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
        $rH_retVal = _pairCommand();
    }
    else {
        die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
    }

    return $rH_retVal;
}

sub _pairCommand {
    my %retVal;
    my $minQuality = $rH_cfg->{'trim.minQuality'};
    my $minLength  = $rH_cfg->{'trim.minLength'};

    my $command             = '';
    my $laneDirectory       = "reads/";
    my $inputFastqPair1Name = $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair1.fastq.gz';
    my $inputFastqPair2Name = $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair2.fastq.gz';
    my $outputFastqPair1Name = $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair1.dup.fastq.gz';
    my $outputFastqPair2Name = $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair2.dup.fastq.gz';
    my $pair1FileDate = -M $outputFastqPair1Name;
    my $pair2FileDate = -M $outputFastqPair2Name;

    my $currentFileDate = $pair2FileDate;
    if ( defined($pair1FileDate) && $pair1FileDate > $pair2FileDate ) {
        $currentFileDate = $pair1FileDate;
    }

    if ( !defined($currentFileDate) || $currentFileDate > -M $rH_laneInfo->{'read1File'} ) {
        $command .= 'module add jdk;';
        $command .= ' ' . $rH_cfg->{'duplicate.javaProgram'} . ' -i ' . $inputFastqPair1Name . ' -i2 ' . $inputFastqPair2Name;
        $command .= ' -k 20 -o 15;';
        $command .= ' mv ' . $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair1.fastq.gz.dup.read1.gz ' . $outputFastqPair1Name . ';';
        $command .= ' mv ' . $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair2.fastq.gz.dup.read1.gz ' . $outputFastqPair2Name . ';';
    }
    $retVal{'command'} = $command;
    $retVal{'pair1'}   = $outputFastqPair1Name;
    $retVal{'pair2'}   = $outputFastqPair2Name;

    return ( \%retVal );

}

sub _singleCommand {
    my %retVal;
    my $minQuality = $rH_cfg->{'trim.minQuality'};
    my $minLength  = $rH_cfg->{'trim.minLength'};

    my $command         = '';
    my $laneDirectory   = "reads/";
    my $inputFastqName  = $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single.fastq.gz';
    my $outputFastqName = $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single.dup.fastq.gz';
    my $currentFileDate = -M $outputFastqName;

    if ( $currentFileDate > -M $rH_laneInfo->{'read1File'} ) {
        $command .= 'module add jdk;';
        $command .= ' ' . $rH_cfg->{'duplicate.javaProgram'} . ' -i ' . $inputFastqName;
        $command .= ' -k 20 -o 15;';
        $command .= ' mv ' . $laneDirectory . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single.fastq.gz.dup.read1.gz ' . $outputFastqName . ';';
    }
    $retVal{'command'} = $command;
    $retVal{'single1'} = $outputFastqName;

    return ( \%retVal );
}

1;

