#!/usr/env/perl

=head1 NAME

I<Trimmomatic>

=head1 SYNOPSIS

Trimmomatic->trim()

=head1 DESCRIPTION

B<Trimmomatic> is a library that trims fastqs

Input = file_name

Output = array


=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package Trimmomatic;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub trim {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;
    my $outputDir   = shift;

    my $rH_retVal;

    if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
        $rH_retVal = singleCommand( $rH_cfg, $sampleName, $rH_laneInfo, $outputDir );
    }
    elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
        $rH_retVal = pairCommand( $rH_cfg, $sampleName, $rH_laneInfo, $outputDir );
    }
    else {
        die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
    }

    return $rH_retVal;
}

sub pairCommand {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;
    my $outputDir   = shift;

    my $minQuality  = LoadConfig::getParam($rH_cfg, 'trim','minQuality');
    my $minLength   = LoadConfig::getParam($rH_cfg, 'trim','minLength');
    my $adapterFile = LoadConfig::getParam($rH_cfg, 'trim','adapterFile');

    my $rawReadDir    = LoadConfig::getParam($rH_cfg, 'trim','rawReadDir');

    my $outputFastqPair1Name = $outputDir .'/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair1.fastq.gz';
    my $outputFastqPair2Name = $outputDir .'/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair2.fastq.gz';
    my $outputFastqSingle1Name = $outputDir .'/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single1.fastq.gz';
    my $outputFastqSingle2Name = $outputDir .'/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single2.fastq.gz';
    my $outputTrimLog = $outputDir .'/' . $sampleName . '.trim.out';
    my $pair1FileDate = -M $outputFastqPair1Name;
    my $pair2FileDate = -M $outputFastqPair2Name;

    my $currentFileDate = $pair2FileDate;
    if ( defined($pair1FileDate) && $pair1FileDate > $pair2FileDate ) {
        $currentFileDate = $pair1FileDate;
    }

    my $command = "";

    # -M gives modified date relative to now. The bigger the older.
    if ( !defined($currentFileDate) || $currentFileDate > -M $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'} ) {
        $command .= 'module load';
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'trim','moduleVersion.java');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'trim','moduleVersion.trimmomatic');
        $command .= ' ; java -cp \$TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE';
        $command .= ' -threads ' . $rH_cfg->{'trim.nbThreads'};
        if ( $rH_laneInfo->{'qualOffset'} eq "64" ) {
            $command .= ' -phred64';
        }
        else {
            $command .= ' -phred33';
        }
        $command .= ' ' . $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'} . ' ' . $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read2File'};
        $command .= ' ' . $outputFastqPair1Name . ' ' . $outputFastqSingle1Name;
        $command .= ' ' . $outputFastqPair2Name . ' ' . $outputFastqSingle2Name;
        if ( $rH_laneInfo->{'qualOffset'} eq "64" ) {
            $command .= ' TOPHRED33';
        }
        $command .= ' ILLUMINACLIP:' . $adapterFile . $rH_cfg->{'trim.clipSettings'};
        if ( $minQuality > 0 ) {
            $command .= ' TRAILING:' . $minQuality;
        }
        $command .= ' MINLEN:' . $minLength;
        $command .= ' > ' . $outputTrimLog;
    }

    my %retVal;
    $retVal{'command'} = $command;
    $retVal{'pair1'}   = $outputFastqPair1Name;
    $retVal{'pair2'}   = $outputFastqPair2Name;
    $retVal{'single1'} = $outputFastqSingle1Name;
    $retVal{'single2'} = $outputFastqSingle2Name;
    return \%retVal;
}

sub singleCommand {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;
    my $outputDir   = shift;

    my $minQuality  = LoadConfig::getParam($rH_cfg, 'trim','minQuality');
    my $minLength   = LoadConfig::getParam($rH_cfg, 'trim','minLength');
    my $adapterFile = LoadConfig::getParam($rH_cfg, 'trim','adapterFile');

    my $rawReadDir    = LoadConfig::getParam($rH_cfg, 'trim','rawReadDir');

    my $outputFastqName = $outputDir . '/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single.fastq.gz';
    my $outputTrimLog = $outputDir . '/' . $sampleName . '.trim.out';
    my $currentFileDate = -M $outputFastqName;

    my $command = "";

    # -M gives modified date relative to now. The bigger the older.
    if ( !defined($currentFileDate) || $currentFileDate > -M $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'} ) {
        $command .= 'module load';
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'trim','moduleVersion.java');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'trim','moduleVersion.trimmomatic');
        $command .= ' ; java -cp \$TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticSE';
        $command .= ' -threads ' . $rH_cfg->{'trim.nbThreads'};
        if ( $rH_laneInfo->{'qualOffset'} eq "64" ) {
            $command .= ' -phred64';
        }
        else {
            $command .= ' -phred33';
        }
        $command .= ' ' . $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'} . ' ' . $outputFastqName;
        $command .= ' ILLUMINACLIP:' . $adapterFile . ':2:30:15';
        if ( $minQuality > 0 ) {
            $command .= ' TRAILING:' . $minQuality;
        }
        $command .= ' MINLEN:' . $minLength;
        $command .= ' > ' . $outputTrimLog;
    }

    my %retVal;
    $retVal{'command'} = $command;
    $retVal{'single1'} = $outputFastqName;
    return \%retVal;
}

1;
