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
        my $defaultDir;
        my $qualOffSet;
        foreach my $sampleName ( keys %{$rHoA0H_sampleSheetInfo} ) {
            my $rAoH_sampleLanes = $rHoA0H_sampleSheetInfo->{$sampleName};
            foreach my $rH_laneInfo (@$rAoH_sampleLanes) {
                if ( $line[0] eq $rH_laneInfo->{'name'} ) {

                    $defaultDir = $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '/';
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
        $sampleInfo{'runType'}                   =  $pairType; 

        #Pair 1
        if ( $pairType eq "PAIRED_END" ) {

            #----------------------------------------------------------------
            # $groupInfo{genome}{-left} =  sample1 t sample2 ...
            # $groupInfo{genome}{-right} =  sample1  sample2 ..
            #-----------------------------------------------------------------
            $groupInfo{'group'}{$line[1] }->{'left'} .= '  ' . $sampleInfo{ $line[0] }->{'sample_name'} . '/' . $defaultDir . $sampleInfo{ $line[0] }->{'sample_name'} . '.t' . $minQuality . 'l' . $minLength .  '.pair1.fastq.gz.dup.gz ';
            $groupInfo{'group'}{ $line[1] }->{'right'} .= '  ' . $sampleInfo{ $line[0] }->{'sample_name'} . '/' . $defaultDir . $sampleInfo{ $line[0] }->{'sample_name'} . '.t' . $minQuality . 'l' . $minLength . '.pair1.fastq.gz.dup.gz ';


        }
        elsif ( $pairType eq "SINGLE_END"){
        	$groupInfo{'group'}{$line[1] }->{'single'} .= ' --single ' .  $sampleInfo{ $line[0] }->{'sample_name'} . '/' . $defaultDir . $sampleInfo{ $line[0] }->{'sample_name'} . '.t' . $minQuality . 'l' . $minLength . '.single.fastq.gz';
        }

    }

    return ( \%sampleInfo, \%groupInfo );

}

1;

