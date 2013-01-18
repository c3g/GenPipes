#!/usr/env/perl

=head1 NAME

I<SplitFile>

=head1 SYNOPSIS

SplitFile::sub(args)

B<Splitting a multiFasta>

SplitFile::splitFasta($fileName, %ref_hash_config, $sample_name, %ref_hash_laneInfo)


B<Splitting Butterfly commands>

SplitFile::splitButterfly(%ref_hash_config, $sample_name, %ref_hash_laneInfo)


=head1 DESCRIPTION

B<SplitFile> is a library used to split files in
a custom way. Each sub is resposible for a type of
file.


=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug
=cut

package SplitFile;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Data::Dumper;
use Config::Simple;

#---------------
# SUB
#-------------

our $rH_cfg;
our $sampleName;
our $rH_laneInfo;

sub splitFasta {
    my $fileName = shift;
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;

    my %retVal;
    my $command        = '';
    my $laneDirectory  = 'alignment/' . $sampleName . "/";
    my $splitDirectory = $laneDirectory . "fasta_split/";

    $command .= ' mkdir -p ' . $splitDirectory . ' ; ';
    $command .= $rH_cfg->{'blast.split'} . ' -c ' . $rH_cfg->{'blast.chunks'};
    $command .= ' -f ' . $laneDirectory . $fileName . ' -o ' . $splitDirectory;

    $retVal{'command'} = $command;
    return ( \%retVal );
}

sub splitButterfly {

    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;

    my %retVal;
    my $command        = '';
    my $laneDirectory  = 'alignment/' . $sampleName . '/chrysalis/';
    my $splitDirectory = $laneDirectory . "butterfly_split";

    $command .= ' mkdir -p ' . $splitDirectory . ' ; ';
    $command .= ' split -a 4 -d -l ' . $rH_cfg->{'trinity.splitLines'} . ' ' . $laneDirectory . 'butterfly_commands.adj '; # split numerically with a 4 digits padding
    $command .= $splitDirectory . '/cmd.';

    $retVal{'command'} = $command;
    return ( \%retVal );

}

1;

