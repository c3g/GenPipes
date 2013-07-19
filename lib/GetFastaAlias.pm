#!/usr/env/perl

=head1 NAME

I<GetFastaAlias>

=head1 SYNOPSIS

GetFastaAlias::sampleInfo(File_path, %Config_hash_ref, %rHoA0H_sampleSheetInfo_ref)

=head1 DESCRIPTION

B<GetFastaAlias> is a library used by De Novo transcriptome assembly,
which retrieves information about samples and genome groups.

It returns two reference hashes (%sampleInfo and %groupInfo);


=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug



=cut

package GetFastaAlias;

# Strict Pragmas
#--------------------------
use Cwd;
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Data::Dumper;
use Config::Simple;

sub sampleInfo {
    my $file                   = shift;
    my $rH_cfg                 = shift;
    my $rHoA0H_sampleSheetInfo = shift;
    my $pairType;
    my %sampleInfo;
    my %groupInfo;

    open( IN, $file ) || die "could not open $file, $!\n";
    while (<IN>) {
        chomp;
        my @line = split /\,/, $_;
        my $defaultDir;               # Directory of the trimmed reads
        my $qualOffSet;
        foreach my $sampleName ( keys %{$rHoA0H_sampleSheetInfo} ) {
            my $rAoH_sampleLanes = $rHoA0H_sampleSheetInfo->{$sampleName};
            foreach my $rH_laneInfo (@$rAoH_sampleLanes) {
                if ( $line[0] eq $rH_laneInfo->{'name'} ) {

                    $defaultDir = getcwd().'/reads/';
                    $qualOffSet = $rH_laneInfo->{'qualOffset'};
                    $pairType   = $rH_laneInfo->{'runType'};

                    if ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
                        $sampleInfo{ $line[0] }->{'pair1'} = $rH_laneInfo->{'read1File'};
                        $sampleInfo{ $line[0] }->{'pair2'} = $rH_laneInfo->{'read2File'};
                    }
                    else {
                        $sampleInfo{ $line[0] }->{'single'} = $rH_laneInfo->{'read1File'};
                    }

                }
            }
        }

        my $minQuality = $rH_cfg->{'trim.minQuality'};
        my $minLength  = $rH_cfg->{'trim.minLength'};

        $sampleInfo{ $line[0] }->{'path'}        = $line[0] . '/' . $defaultDir;
        $sampleInfo{ $line[0] }->{'group_name'}  = $line[1];
        $sampleInfo{ $line[0] }->{'sample_name'} = $line[0];
        $sampleInfo{'defaultDir'}                = $defaultDir;
        $sampleInfo{'runType'}                   = $pairType;

        if(!defined($groupInfo{ $line[1] })){
          $groupInfo{ $line[1] } = {};
        }

        #Pair 1
        if ( $pairType eq "PAIRED_END" ) {

            #----------------------------------------------------------------
            # $groupInfo{genome}->{-left} =  sample1 t sample2 ...
            # $groupInfo{genome}->{-right} =  sample1  sample2 ..
            #-----------------------------------------------------------------
            $groupInfo{ $line[1] }->{'left'} .= ' '  . $sampleInfo{'defaultDir'} . $line[0] . '/' . $sampleInfo{ $line[0] }->{'sample_name'} . '.pair1.dup.fastq.gz ';
            $groupInfo{ $line[1] }->{'right'} .= ' ' . $sampleInfo{'defaultDir'} . $line[0] . '/' . $sampleInfo{ $line[0] }->{'sample_name'} . '.pair2.dup.fastq.gz ';

        }
        elsif ( $pairType eq "SINGLE_END" ) {
            $groupInfo{ $line[1] }->{'single'} .= ' --single ' . $defaultDir . $sampleInfo{ $line[0] }->{'sample_name'} . '.single.dup.fastq.gz';
        }

    }

    return ( \%sampleInfo, \%groupInfo );

}

1;

