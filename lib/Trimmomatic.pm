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
use PipelineUtils;
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

    my $ro_job;

    if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
        $ro_job = singleCommand( $rH_cfg, $sampleName, $rH_laneInfo, $outputDir );
    }
    elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
        $ro_job = pairCommand( $rH_cfg, $sampleName, $rH_laneInfo, $outputDir );
    }
    else {
        die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
    }

    return $ro_job;
}

sub pairCommand {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;
    my $outputDir   = shift;

    my $minQuality  = LoadConfig::getParam($rH_cfg, 'trim','minQuality');
    my $minLength   = LoadConfig::getParam($rH_cfg, 'trim','minLength');
    my $adapterFile = LoadConfig::getParam($rH_cfg, 'trim','adapterFile');
    my $headcrop    = LoadConfig::getParam($rH_cfg, 'trim','headcrop');

    my $rawReadDir    = LoadConfig::getParam($rH_cfg, 'trim','rawReadDir');

    my $outputFastqPair1Name = $outputDir .'/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair1.fastq.gz';
    my $outputFastqPair2Name = $outputDir .'/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.pair2.fastq.gz';
    my $outputFastqSingle1Name = $outputDir .'/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single1.fastq.gz';
    my $outputFastqSingle2Name = $outputDir .'/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single2.fastq.gz';
    my $outputTrimLog = $outputDir .'/' . $sampleName . '.trim.out';
    my $outputTrimStats = $outputDir .'/' . $sampleName . '.trim.stats.csv';

    my $inputFastqPair1Name = $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};
    my $inputFastqPair2Name = $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read2File'};

    my $up2date = PipelineUtils::testInputOutputs([$inputFastqPair1Name, $inputFastqPair2Name], [$outputFastqPair1Name,$outputFastqPair2Name,$outputFastqSingle1Name,$outputFastqSingle2Name,$outputTrimLog,$outputTrimStats]);

    my $ro_job = new Job(!defined($up2date));

    # -M gives modified date relative to now. The bigger the older.
    if (!$ro_job->isUp2Date()) {
        my $command = "";
        $command .= 'module load';
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'trim','moduleVersion.java');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'trim','moduleVersion.trimmomatic');
        $command .= ' ; java -XX:ParallelGCThreads=1 -Xmx2G -cp \$TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE';
        $command .= ' -threads ' . $rH_cfg->{'trim.nbThreads'};
        if ( $rH_laneInfo->{'qualOffset'} eq "64" ) {
            $command .= ' -phred64';
        }
        else {
            $command .= ' -phred33';
        }
        $command .= ' ' . $inputFastqPair1Name;
        $command .= ' ' . $inputFastqPair2Name;
        $command .= ' ' . $outputFastqPair1Name . ' ' . $outputFastqSingle1Name;
        $command .= ' ' . $outputFastqPair2Name . ' ' . $outputFastqSingle2Name;
        if ( $rH_laneInfo->{'qualOffset'} eq "64" ) {
            $command .= ' TOPHRED33';
        }
        if(defined($headcrop) && length($headcrop) > 0 && $headcrop > 0) {
          $command .= ' HEADCROP:' . $headcrop;
        }
        $command .= ' ILLUMINACLIP:' . $adapterFile . $rH_cfg->{'trim.clipSettings'};
        if ( $minQuality > 0 ) {
            $command .= ' TRAILING:' . $minQuality;
        }
        $command .= ' MINLEN:' . $minLength;
        $command .= ' 2> ' . $outputTrimLog;
        $command .= ' &&';
        $command .= ' grep \"^Input Read\" '.$outputTrimLog.'| sed \'s/Input Read Pairs: \\([0-9]\\+\\).*Both Surviving: \\([0-9]\\+\\).*Forward Only Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\3/g\' | tr \'#\' \'\n\' > '.$outputTrimStats;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }

    return $ro_job;
}

sub singleCommand {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;
    my $outputDir   = shift;

    my $minQuality  = LoadConfig::getParam($rH_cfg, 'trim','minQuality');
    my $minLength   = LoadConfig::getParam($rH_cfg, 'trim','minLength');
    my $adapterFile = LoadConfig::getParam($rH_cfg, 'trim','adapterFile');
    my $headcrop    = LoadConfig::getParam($rH_cfg, 'trim','headcrop');

    my $rawReadDir    = LoadConfig::getParam($rH_cfg, 'trim','rawReadDir');

    my $outputFastqName = $outputDir . '/' . $sampleName . '.t' . $minQuality . 'l' . $minLength . '.single.fastq.gz';
    my $outputTrimLog = $outputDir . '/' . $sampleName . '.trim.out';
    my $outputTrimStats = $outputDir .'/' . $sampleName . '.trim.stats.csv';
    my $inputFastqName = $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};

    my $up2date = PipelineUtils::testInputOutputs([$inputFastqName], [$outputFastqName,$outputTrimLog,$outputTrimStats]);

    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command = "";
        $command .= 'module load';
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'trim','moduleVersion.java');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'trim','moduleVersion.trimmomatic');
        $command .= ' ; java -XX:ParallelGCThreads=1 -Xmx2G -cp \$TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticSE';
        $command .= ' -threads ' . $rH_cfg->{'trim.nbThreads'};
        if ( $rH_laneInfo->{'qualOffset'} eq "64" ) {
            $command .= ' -phred64';
        }
        else {
            $command .= ' -phred33';
        }
        $command .= ' ' . $inputFastqName;
        $command .= ' ' . $outputFastqName;
        if ( $rH_laneInfo->{'qualOffset'} eq "64" ) {
            $command .= ' TOPHRED33';
        }
        if(defined($headcrop) && length($headcrop) > 0 && $headcrop > 0) {
          $command .= ' HEADCROP:' . $headcrop;
        }
        $command .= ' ILLUMINACLIP:' . $adapterFile . $rH_cfg->{'trim.clipSettings'};
        if ( $minQuality > 0 ) {
            $command .= ' TRAILING:' . $minQuality;
        }
        $command .= ' MINLEN:' . $minLength;
        $command .= ' 2> ' . $outputTrimLog;
        $command .= ' &&';
        $command .= ' grep \"^Input Read\" '.$outputTrimLog.'| sed \'s/Input Reads: \\([0-9]\\+\\).*Surviving: \\([0-9]\\+\\).*/Raw Fragments,\\1#Fragment Surviving,\\2#Single Surviving,\\2/g\' | tr \'#\' \'\n\' > '.$outputTrimStats;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }

    return $ro_job;
}

1;
