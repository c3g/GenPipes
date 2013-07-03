#!/usr/env/perl

=head1 NAME

I<Trinity>

=head1 SYNOPSIS

Trinity::sub(args)

B<Trinity::chrysalis>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $pair1, $pair2). In the case where there are multiple samples
the $pair1 and $pair2 will be  string with all the  samples (reads/sample1_left   reads/sample2_left).
 
B<Trinity::butterfly>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $fileButterflyComand)

B<Trinity::concatFastaCreateGtf>(%ref_hash_config, $sample_name, %ref_hash_laneInfo)


All subroutines return a ref_hash with the command line



=head1 DESCRIPTION

B<Trinity> is a library to use the
transcriptome assembly package, Trinity.

The lib has three subroutines: B<chrysalis, butterfly and concatFastaCreateGtf>.

- chrysalis and butterfly are steps of the package trinity.

- concatFastaCreateGtf is used to concatenate the fasta files and 
to generate de gtf file.

Each sub must be called independently with their set of variables.

=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug



=cut

package Trinity;

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
our $pair1;
our $pair2;
our $fileButterflyComand;

sub chrysalis {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $pair1       = shift;    # For single command the left will receive the file.
    $pair2       = shift;

    my $rH_retVal;
    if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
        $rH_retVal = _chrysalisSingleCommand();
    }
    elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
        $rH_retVal = _chrysalisPairCommand();
    }
    else {
        die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
    }

}

sub butterfly {
    $rH_cfg              = shift;
    $sampleName          = shift;
    $rH_laneInfo         = shift;
    $fileButterflyComand = shift;

    my $laneDirectory = "assembly/" . $sampleName . "/chrysalis/";
    my $command       = ' ';
    my %retVal;

    $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
    $command .= ' ' . $rH_cfg->{'butterfly.parallel'} . ' -f ' . $laneDirectory . 'butterfly_split/' . $fileButterflyComand;
    $command .= ' -n ' . $rH_cfg->{'butterfly.nbThreads'} . ' ';

    $retVal{'command'} = $command;
    return ( \%retVal );

}

sub concatFastaCreateGtf {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;

    my $command = '';
    my %retVal;
    my $laneDirectory = "assembly/" . $sampleName . "/";

    $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
    $command .= ' find ' . $laneDirectory . 'chrysalis';
    $command .= ' -name "*allProbPaths.fasta" -exec cat {} + >' . $laneDirectory . 'Trinity.fasta ;';
    $command .= ' sh ' . $rH_cfg->{'trinity.createGtf'} . ' ' . $laneDirectory . 'Trinity.fasta';
    $command .= ' ' . $laneDirectory . $sampleName . '.gtf ;';
    $command .= ' awk \'{print \$1} \' ' . $laneDirectory . 'Trinity.fasta ';
    $command .= ' >' . $laneDirectory . 'Trinity.2.fasta';

    $retVal{'command'} = $command;
    return ( \%retVal );

}

sub _chrysalisPairCommand {

    my $command = '';
    my %retVal;

    my $laneDirectory = 'assembly/' . $sampleName . '/';
    $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
    $command .= ' Trinity.pl --seqType fq --JM 100G';
    $command .= ' --left' . ' \" ' . $pair1 . ' \" ' . '--right' . ' \" ' . $pair2 . ' \" ';
    $command .= ' --CPU ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'nbThreads');
    $command .= ' --output ' . $laneDirectory;
    $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

    $retVal{'command'} = $command;
    return ( \%retVal );

}

sub _chrysalisSingleCommand {

    my $command = '';
    my %retVal;
    
    my $laneDirectory = 'assembly/' . $sampleName . '/';
    $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
    $command .= ' Trinity.pl --seqType fq --JM 100G';
    $command .= ' ' . $pair1;
    $command .= ' --CPU ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'nbThreads');
    $command .= ' --output ' . $laneDirectory;
    $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

    $retVal{'command'} = $command;
    return ( \%retVal );

}

1;

