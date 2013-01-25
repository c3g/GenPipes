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
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

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
    
    my $laneDirectory = 'DGE/';
	my $command = '';
	
	
	$command .= ' module add gcc/4.7.0 ; module add R ;';
	$command .= ' mkdir -p DGE/' . $group . ';';
	$command .= ' Rscript ' . $rH_cfg->{'diffExpress.edger'} . ' -d ' . $rH_cfg->{'diffExpress.designFile'} ;
	$command .= ' -c ' .  $laneDirectory . $group . '/matrix.csv ';
	$command .= ' -o ' . $laneDirectory . $group . '/ ;' ;
	$command .= ' Rscript ' . $rH_cfg->{'diffExpress.deseq'} . ' -d ' . $rH_cfg->{'diffExpress.designFile'} ;
	$command .= ' -c ' .  $laneDirectory . $group . '/matrix.csv ';
	$command .= ' -o ' . $laneDirectory . $group . '/ ;' ;
	
	$retVal{'command'} = $command;
	return ( \%retVal );
		
	
}

sub edgerPortable {
	my $rH_cfg        = shift;
	my $designFile    = shift;
	my $countMatrix   = shift;
	my $outputDir     = shift;
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'diffExpress','toolsModule') .' ' .LoadConfig::getParam($rH_cfg, 'diffExpress','cranRModule') .' ;';
	$command .= ' Rscript $R_TOOLS/edger.R -d' .$designFile;
	$command .= ' -c' .$countMatrix;
	$command .= ' -o' .$outputDir;

	return $command;
}

sub deseq {
	my $rH_cfg        = shift;
	my $designFile    = shift;
	my $countMatrix   = shift;
	my $outputDir     = shift;
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'diffExpress','toolsModule') .' ' .LoadConfig::getParam($rH_cfg, 'diffExpress','cranRModule') .' ;';
	$command .= ' Rscript $R_TOOLS/deseq.R -d' .$designFile;
	$command .= ' -c' .$countMatrix;
	$command .= ' -o' .$outputDir;

	return $command;
}

sub matrix {
	my $rH_cfg        = shift;
	
}

1;

