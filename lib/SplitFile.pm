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
sub splitFasta {
    my $fileName = shift;
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;

    my %retVal;
    my $laneDirectory  = 'assembly/' . $sampleName . "/";
    my $splitDirectory = $laneDirectory . "fasta_split/";

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$laneDirectory . $fileName],undef);

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= ' mkdir -p ' . $splitDirectory . ' &&';
      $command .= ' module load '.LoadConfig::getParam($rH_cfg, 'blast', 'moduleVersion.exonerate').' &&';
      $command .= ' fastasplit' . ' -c ' . $rH_cfg->{'blast.chunks'};
      $command .= ' -f ' . $laneDirectory . $fileName . ' -o ' . $splitDirectory;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub splitButterfly {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;

    my %retVal;
    my $laneDirectory  = 'assembly/' . $sampleName . '/chrysalis/';
    my $splitDirectory = $laneDirectory . "butterfly_split";

    my $ro_job = new Job();
    $ro_job->testInputOutputs(undef,undef);

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= ' mkdir -p ' . $splitDirectory . ' && ';
      $command .=  LoadConfig::getParam($rH_cfg, 'butterfly', 'split') . ' ' . LoadConfig::getParam($rH_cfg, 'butterfly', 'chunks') . ' ' . $laneDirectory . 'butterfly_commands.adj '; # split numerically with a 4 digits padding
      $command .= $splitDirectory . ' ';

      $ro_job->addCommand($command);
    }
    return $ro_job;

}

1;

