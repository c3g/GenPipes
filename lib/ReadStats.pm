#!/usr/env/perl

=head1 NAME

I<ReadStats>

=head1 SYNOPSIS

ReadStats::sub(args)

B<ReadStats::stats>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $read, $sortedBam, $group) 

B<ReadStats::concatStats>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group) 

B<ReadStats::contigStats>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group)

All subroutines return a ref_hash with the command line


=head1 DESCRIPTION

B<readStats> is a library to generate basic statistics of each sample/contig.

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
our $read;
our $sortedBam;
our $group;

sub stats {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $read        = shift;
    $sortedBam   = shift;
    $group       = shift;

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

sub _pairCommand {
    my %retVal;
    my $command = '';

    $group = ( !defined $group ) ? $sampleName : $group;

    my $laneDirectory = "stats/" . $group . "/";
    my $minQuality    = $rH_cfg->{'trim.minQuality'};
    my $minLength     = $rH_cfg->{'trim.minLength'};

    $command .= ' module add ' . LoadConfig::getParam($rH_cfg, 'readstats','moduleVersion.samtools') . ' ;';
    $command .= ' mkdir -p ' . $laneDirectory . ';';
    $command .= ' zcat ' . $rH_cfg->{'default.rawReadDir'} . '/' . $sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};
    $command .= ' | wc -l | awk \'{print \$1}\' | echo \$[\$(</dev/stdin)/4]';
    $command .= ' >' . $laneDirectory . $sampleName . '.readstats.raw.csv ; ';
    $command .= ' zcat ' . $read;
    $command .= '| wc -l | awk \'{print \$1}\' | echo \$[\$(</dev/stdin)/4]';
    $command .= ' >' . $laneDirectory . $sampleName . '.readstats.filtered.csv ;';
    $command .= ' samtools view -F4 -F256 ' . $sortedBam . ' | ';
    $command .= ' grep -v \"^[@]\" | awk \'BEGIN={ si=0;pa=0 } {if (\$3 !="*" ) ';
    $command .= ' if {(\$7=="*" ) {pa=pa+1} else{ si=si+1 }}} END { print \$si+(pa/1) }\'';
    $command .= ' >' . $laneDirectory . $sampleName . '.readstats.aligned.csv';

    $retVal{'command'} = $command;
    return ( \%retVal );

}

sub _singleCommand {
    my %retVal;
    my $command = '';

    $group = ( !defined $group ) ? $sampleName : $group;

    my $laneDirectory = "stats/" . $group . "/";
    my $minQuality    = $rH_cfg->{'trim.minQuality'};
    my $minLength     = $rH_cfg->{'trim.minLength'};
    
    $command .= ' module add ' . LoadConfig::getParam($rH_cfg, 'readstats','moduleVersion.samtools') . ' ;';
    $command .= ' mkdir -p ' . $laneDirectory . ';';
    $command .= ' zcat ' . $rH_cfg->{'default.rawReadDir'} . '/' . $rH_laneInfo->{'read1File'};
    $command .= ' | wc -l | awk \'{print$1}\' | echo $[$(</dev/stdin)/4]';
    $command .= ' >' . $laneDirectory . $sampleName . '.readstats.raw.csv ; ';
    $command .= ' zcat ' . $read;
    $command .= '| wc -l | awk \'{print$1}\' | echo $[$(</dev/stdin)/4]';
    $command .= ' >' . $laneDirectory . $sampleName . '.readstats.filtered.csv ;';
    $command .= ' samtools view -F4 -F256 ' . $sortedBam . ' | ';
    $command .= ' grep -v \"^[@]\" | awk \'BEGIN={si=0} {if ($3 !="*" ) ';
    $command .= ' si=si+1} END {print $si}\'';
    $command .= ' >' . $laneDirectory . $sampleName . '.readstats.aligned.csv';

    $retVal{'command'} = $command;
    return ( \%retVal );
}

sub concatStats {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;

    $group = ( !defined $group ) ? $sampleName : $group;
    my $laneDirectory = "stats/" . $group . "/";

    my %retVal;
    my $command = '';

    $command .= ' echo ' . $sampleName . ' > ' . $laneDirectory . $sampleName . '_temp.txt ;';
    $command .= ' paste -d "," ' . $laneDirectory . $sampleName . '_temp.txt ' . $laneDirectory . $sampleName . '.readstats.raw.csv  ';
    $command .= $laneDirectory . $sampleName . '.readstats.filtered.csv ' . $laneDirectory . $sampleName . '.readstats.aligned.csv';
    $command .= ' >>' . $laneDirectory . 'allSamples.' . $group . '.readstats.csv ;';

    $retVal{'command'} = $command;
    return ( \%retVal );
}

sub contigStats {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;

    $group = ( !defined $group ) ? $sampleName : $group;
    my $laneDirectory = "stats/" . $group . "/";

    my %retVal;
    my $command = '';

    $command .= $rH_cfg->{'readStats.fastaLength'} . ' assembly/Trinity.fasta';
    $command .= ' | sort -k1,1nr >' . $laneDirectory . $group . '.contig_length_sorted.txt ;';
    $command .= ' sh ' . $rH_cfg->{'readStats.contigStat'} . ' ' . $laneDirectory . $group . '.contig_length_sorted.txt';
    $command .= ' > ' . $laneDirectory . 'contigs.' . $group . '.stats.txt ;';

    $retVal{'command'} = $command;
    return ( \%retVal );
}

1;

