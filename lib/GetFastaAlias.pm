#!/usr/env/perl

=head1 NAME

I<GetFastaAlias>

=head1 SYNOPSIS

GetFastaAlias::sampleInfo(File_path, %Config_hash, $baseOutput_Dir)

=head1 DESCRIPTION

B<GetFastaAlias> is a library used by De Novo transcriptome assembly,
which retrieves information about samples and genome groups.

It returns two reference hashes (%sampleInfo and %groupInfo);


=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<File::Basename> Path parsing

B<Cwd> Path parsing

=cut

package GetFastaAlias;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Data::Dumper;
use Config::Simple;
use File::Basename;
use Cwd 'abs_path';

sub sampleInfo {
    my $file        = shift;
    my $rH_cfg      = shift;
    my $outDir      = shift;
    my $pairCounter = 1;
    my %sampleInfo;
    my %groupInfo;

    open( IN, $file ) || die "could not open $file, $!\n";
    while (<IN>) {
        chomp;
        my @line = split /\,/, $_;

        $sampleInfo{ $line[1] }->{"path"}        = dirname( $line[0] );
        $sampleInfo{ $line[1] }->{"group_name"}  = $line[2];
        $sampleInfo{ $line[1] }->{"sample_name"} = $line[1];

        #Pair 1
        if ( $pairCounter == 1 ) {
            $sampleInfo{ $line[1] }->{"pair1"} = basename( $line[0] );

            #----------------------------------------------------------------
            # $groupInfo{genome}{-left} = --left sample1 --left sample2 ...
            #-----------------------------------------------------------------
            $groupInfo{ $line[2] }->{"left"} .= " --left " . "$outDir/reads/$sampleInfo{$line[1]}->{'sample_name'}_t$rH_cfg->{'trim.minQuality'}l$rH_cfg->{'trim.minLength'}.phred33.pair1.fastq.gz";
            $pairCounter = 2;
        }

        # Pair 2
        elsif ( $pairCounter == 2 )
        {
            $sampleInfo{ $line[1] }->{"pair2"} = basename( $line[0] );

            #----------------------------------------------------------------
            # $groupInfo{genome}{-right} = --right sample1 --right sample2 ...
            #-----------------------------------------------------------------
            $groupInfo{ $line[2] }->{"right"} .= " --right " . "$outDir/reads/$sampleInfo{$line[1]}->{'sample_name'}_t$rH_cfg->{'trim.minQuality'}l$rH_cfg->{'trim.minLength'}.phred33.pair2.fastq.gz";
            $pairCounter = 1;

        }

        # Single end
        else {
            $sampleInfo{ $line[1] }->{"single"} = basename( $line[0] );
        }

    }

    return ( \%sampleInfo, \%groupInfo );

}

1;

