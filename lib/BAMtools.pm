#!/usr/env/perl

=head1 NAME

I<SVTools>

=head1 SYNOPSIS

Picard->merge()

=head1 DESCRIPTION

B<SVTools> is a library to analyse BAMs for Structural Variants

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package BAMtools;

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
sub countBins {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $tumorBam    = shift;
    my $window      = shift;
    my $normType    = shift;
    my $outputFile  = shift;
    my $normalBam   = shift; # can be undef for non paired run

    my $up2date;
    if(defined($normalBam)) {
        $up2date = PipelineUtils::testInputOutputs([$tumorBam, $normalBam], [$outputFile]);
    }
    else {
        $up2date = PipelineUtils::testInputOutputs([$tumorBam], [$outputFile]);
    }

    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'countBins', 'moduleVersion.java').' ' .LoadConfig::getParam($rH_cfg, 'countBins', 'moduleVersion.bamtools').' ; ';
        $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'countBins', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'countBins', 'extraJavaFlags');
        $command .= ' -Xmx1500M -jar \${BAMTOOLS_JAR} bincounter';
        $command .= ' --norm '.$normType;
        $command .= ' --minMapQ '.LoadConfig::getParam($rH_cfg, 'countBins', 'minMapQ');
        $command .= ' --bam '.$tumorBam;
        if(defined($normalBam)) {
            $command .= ' --refbam '.$normalBam;
        }
        $command .= ' --window '.$window;
        $command .= ' > '.$outputFile;
        $command .= ' '.$up2date;

        $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub deleteDuplicates {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $pair1           = shift;
    my $pair2           = shift;
    my $single          = shift;
    my $optOutputPrefix = shift;


    my $ro_job;
    if (defined($pair1) && defined($pair2)) {
        $ro_job = deletePairedDuplicates($rH_cfg, $sampleName, $pair1, $pair2, $optOutputPrefix);
    }
    elsif (defined($single)) {
        $ro_job = deleteSingleDuplicates($rH_cfg, $sampleName, $single, $optOutputPrefix);
    }
    else {
        die "Unknown runType. \n";
    }

    return $ro_job;
}

sub deletePairedDuplicates {
    my $rH_cfg = shift;
    my $sampleName = shift;
    my $pair1           = shift;
    my $pair2           = shift;
    my $outputPrefix    = shift;
    my %retVal;
	
    my $command             = '';
    my $outputFastqPair1Name = $outputPrefix . '.pair1.dup.fastq.gz';
    my $outputFastqPair2Name = $outputPrefix . '.pair2.dup.fastq.gz';

    my $up2date = PipelineUtils::testInputOutputs([$pair1, $pair2], [$outputFastqPair1Name,$outputFastqPair2Name]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        $command .= 'module add '. LoadConfig::getParam($rH_cfg, 'default','moduleVersion.java') ;
        $command .= ' ' . LoadConfig::getParam($rH_cfg, 'duplicate','moduleVersion.bamtools') . ';' ;
        $command .= ' java '.LoadConfig::getParam($rH_cfg, 'duplicate', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'duplicate', 'dupRam').' -jar  \$BAMTOOLS_JAR filterdups' .  ' --read1 ' . $pair1 . ' --read2 ' . $pair2;
        $command .= ' -k 20 -o 15';
        $command .= ' && mv ' . $pair1 . '.dup.read1.gz ' . $outputFastqPair1Name;
        $command .= ' && mv ' . $pair2 . '.dup.read2.gz ' . $outputFastqPair2Name;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }

    return $ro_job;

}

sub deleteSingleDuplicates {
    my $rH_cfg = shift;
    my $sampleName = shift;
    my $single          = shift;
    my $outputPrefix    = shift;
    my %retVal;
	
    my $outputFastqName = $outputPrefix . '.single.dup.fastq.gz';

    my $up2date = PipelineUtils::testInputOutputs([$single], [$outputFastqName]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'module add '. LoadConfig::getParam($rH_cfg, 'default','moduleVersion.java') ;
        $command .= ' ' . LoadConfig::getParam($rH_cfg, 'duplicate','moduleVersion.bamtools') . ';' ;
        $command .= ' java '.LoadConfig::getParam($rH_cfg, 'duplicate', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'duplicate', 'dupRam').' -jar \$BAMTOOLS_JAR filterdups' . ' --read1 ' . $single;
        $command .= ' -k 20 -o 15';
        $command .= ' && mv ' . $single . '.dup.read1.gz ' . $outputFastqName;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }

    return $ro_job;
}

1;

