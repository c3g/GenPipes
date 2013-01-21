#!/usr/env/perl

=head1 NAME

I<DiffExpression>

=head1 SYNOPSIS

DiffExpression:sub(args)

B<DiffExpression:stats>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group) 

All subroutines return a ref_hash with the command line


=head1 DESCRIPTION

B<DiffExpression> is a library to create differencial expression analysis.

=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug


=cut

package DiffExpression;

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

sub edger{
	$rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
	$group         = shift;
    my %retVal;
    
    my $laneDirectory = 'DGE/' . $group ."/";
	my $command = '';
	
	
	$command .= ' module add gcc/4.7.0 ; module add R ;';
	$command .= ' Rscript ' . $rH_cfg->{'diffExpress.edger'} . ' -d ' . $rH_cfg->{'diffExpress.designFile'} ;
	$command .= ' -c ' .  $laneDirectory . 'matrix.csv ';
	$command .= ' -o ' . $laneDirectory . ' ;' ;
	$command .= ' Rscript ' . $rH_cfg->{'diffExpress.deseq'} . ' -d ' . $rH_cfg->{'diffExpress.designFile'} ;
	$command .= ' -c ' .  $laneDirectory .'matrix.csv ';
	$command .= ' -o ' . $laneDirectory . ' ;' ;
	
	$retVal{'command'} = $command;
	return ( \%retVal );
		
	
}

1;

