#!/usr/env/perl

=head1 NAME

I<ReadStats>

=head1 SYNOPSIS

ReadStats::sub(args)

B<ReadStats::stats>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $readFile, $group) 

B<ReadStats::concatStats>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $readFile, $group) 

All subroutines return a ref_hash with the command line


=head1 DESCRIPTION

B<readStats> is a library to generate basic statistics of each sample.

=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug


=cut

package ReadStats;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Data::Dumper;
use Config::Simple;

our $rH_cfg;
our $sampleName;
our $rH_laneInfo;
our $readFile;
our $group; 

sub stats{
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
	$readFile      = shift;
	$group         = shift;
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
	
	
sub concatStats{
	$rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
	
	my %retVal;
	my $command = '';
	
		
	$command .= ' MISSING COMMAND CONCAT STATS;';
	$retVal{'command'} = $command;
	return ( \%retVal );
}

sub _pairCommand{
	my %retVal;
	my $command = '';
	$group = (!defined $group) ? $sampleName : $group;
	my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/stats" . $group . "/";
	
	$command .= 'mkdir -p' . $laneDirectory . ';';
	$command .= ' MISSING COMMAND ;';
	$retVal{'command'} = $command;
	return ( \%retVal );
	
}

sub _singleCommand{
	my %retVal;
	my $command = '';
	$group = (!defined $group) ? $sampleName : $group;
	my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/stats" . $group . "/";
	
	$command .= 'mkdir -p' . $laneDirectory . ';';
	$command .= ' MISSING SINGLE COMMAND ;';
	$retVal{'command'} = $command;
	return ( \%retVal );
	
}

	
1;










