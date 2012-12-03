#!/usr/env/perl

=head1 NAME

I<Trinity>

=head1 SYNOPSIS

Trinity::sub(args)

B<Trinity::chrysalis>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $pair1, $pair2) if looping per sample

B<Trinity::chrysalis>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, 'undef', 'undef', %ref_hash_groupInfo) if looping per group.
%ref_hash_groupInfo is the ref hash outputed by the package GetFastaAlias (this is probably only useful for transcriptome assembly)
 
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
our $rH_groupInfo;    # used only if assembling multiple samples into one transcriptome (per group)
our $fileButterflyComand;

sub chrysalis {
    $rH_cfg       = shift;
    $sampleName   = shift;
    $rH_laneInfo  = shift;
    $pair1        = shift;
    $pair2        = shift;
    $rH_groupInfo = shift;

    print $sampleName, "SAMPLE \n\n";
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

    my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    my $command       = ' ';
    my %retVal;

    $command .= 'module add jdk64/6u35; ';
    $command .= ' module add mugqic/trinity/2012-06-18 ;';
    $command .= ' ' . $rH_cfg->{'trinity.parallel'} . ' -f ' . $laneDirectory . 'butterfly_split/' . $fileButterflyComand;
    $command .= ' -n ' . $rH_cfg->{'trinity.nbThreads'};

    $retVal{'command'} = $command;
    return ( \%retVal );

}

sub concatFastaCreateGtf {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;

    my $command = '';
    my %retVal;
    my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";

    $command .= 'module add jdk64/6u35; ';
    $command .= ' module add mugqic/trinity/2012-06-18 ;';
    $command .= ' find ' . $laneDirectory . 'chrysalis';
    $command .= ' -name "*allProbPaths.fasta" -exec cat {} + >' . $laneDirectory . 'Trinity.fasta ;';
    $command .= ' sh ' . $rH_cfg->{'trinity.createGtf'} . ' ' . $laneDirectory . 'Trinity.fasta';
    $command .= ' ' . $laneDirectory . $sampleName . '.gtf ;';
    $command .= ' awk \'{print $1}\' ' . $laneDirectory . 'Trinity.fasta ';
    $command .= ' >' . $laneDirectory . 'Trinity.2.fasta';

    $retVal{'command'} = $command;
    return ( \%retVal );

}

sub _chrysalisPairCommand {

    my $command = '';
    my %retVal;

    if ( exists $rH_groupInfo->{'left'} ) {
        my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'defaultDir'};
        $command .= 'module add jdk64/6u35; ';
        $command .= ' module add mugqic/trinity/2012-06-18 ;';
        $command .= ' Trinity.pl --seqType fq --JM 100G';
        $command .= ' ' . $rH_groupInfo->{'left'} . ' ' . $rH_groupInfo->{'right'};
        $command .= ' --CPU 22 --output ' . $laneDirectory;
        $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

    }

    else {
        my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
        $command .= 'module add jdk64/6u35; ';
        $command .= ' module add mugqic/trinity/2012-06-18 ;';
        $command .= ' Trinity.pl --seqType fq --JM 100G';
        $command .= ' --left ' . $pair1 . ' --right ' . $pair2;
        $command .= ' --CPU 22 --output ' . $laneDirectory;
        $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

    }
    $retVal{'command'} = $command;
    return ( \%retVal );

}



sub _chrysalisSingleCommand {

    my $command = '';
    my %retVal;

    if ( exists $rH_groupInfo->{'single'} ) {

        my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
        $command .= 'module add jdk64/6u35; ';
        $command .= ' module add mugqic/trinity/2012-06-18 ;';
        $command .= ' Trinity.pl --seqType fq --JM 100G';
        $command .= ' ' . $rH_groupInfo->{'single'};
        $command .= ' --CPU 22 --output ' . $laneDirectory;
        $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';
    }

    else {
        my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
        $command .= 'module add jdk64/6u35; ';
        $command .= ' module add mugqic/trinity/2012-06-18 ;';
        $command .= ' Trinity.pl --seqType fq --JM 100G';
        $command .= ' --single ' . $pair1;
        $command .= ' --CPU 22 --output ' . $laneDirectory;
        $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

    }

    $retVal{'command'} = $command;
    return ( \%retVal );

}

1;

