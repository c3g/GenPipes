#!/usr/env/perl

=head1 NAME

I<HtseqCount>

=head1 SYNOPSIS

HtseqCount:sub(args)

B<HtseqCount::readCount>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group) 

B<HtseqCount::sortRead>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group)

B<HtseqCount::matrixMake>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group)

All subroutines return a ref_hash with the command line



=head1 DESCRIPTION

B<HtseqCount> is a library to generate basic statistics on raw read count.

=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk
Mathieu Bourgey mbourgey@genomequebec.com

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug


=cut

package HtseqCount;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Data::Dumper;
use Config::Simple;
use SAMtools

our $rH_cfg;
our $sampleName;
our $rH_laneInfo;
our $readFile;
our $group; 

sub matrixMake{
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



sub _singleCommand{
	next;
}


sub _pairCommand{
	next;
}

sub readCount {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $group       = shift;
        $group = ( !defined $group ) ? $sampleName : $group;
        
    my %retVal;
    my $command       = '';
    my $laneDirectory = "read_count/" . $group . "/";

    $command .= ' module add mugqic/samtools/0.1.6; ';
    $command .= ' samtools view ' . $laneDirectory . $sampleName . '.QueryName.bam | ';
    $command .= ' htseq-count - ' . 'alignment/' . $group . '/' . $group . '.gtf ';
    $command .= ' -s no >' . $laneDirectory . $sampleName . '.readcount.cvs';

    $retVal{'command'} = $command;
    return ( \%retVal );
}

sub readCountPortable{
	my $rH_cfg      = shift;
        my $inputBam  = shift;
	my $inputGtf  = shift;
	my $outputFile  = shift;
	my $strandInfo = shift;

	if (!(defined $strandInfo)) {
		$strandInfo='no'
	}
	
	my $command = ;
	my $latestBam = -M $inputBam;
	my $output1 = -M $outputFile;
	if(!defined($latestBam) || !defined($output1) || $latestBam < $output1) {
		$command .= ' module load ' .LoadConfig::getParam($rH_cfg, 'htseq','moduleVersion.python') .' ' .LoadConfig::getParam($rH_cfg, 'htseq','moduleVersion.htseq') .' ; ';
		$command .= ' ; ' .SAMtools::viewFilter($rH_cfg, $inputBam) ;
		$command .= ' | htseq-count - ' .  $inputGtf ;
		$command .= ' -s ' .$strandInfo;
		$command .= ' >' . $outputFile ;
	}
	return $command;
}


sub sortRead{
	$rH_cfg      = shift;
	$sampleName  = shift;
	$rH_laneInfo = shift;
	$group         = shift;
	
	my %retVal;
	my $command = '';
	my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/stats" . $group . "/";
	
	$command .= ' mkdir -p ' . $laneDirectory . 'read_counts/ ;' ;	
	$command .= ' mkdir -p ' . $laneDirectory . 'DGE/ ;' ;	
	$command .= ' module add jdk64 ; module add mugqic/picard/1.64/ ; ';
	$command .= ' java -Xmx30g -jar SortSam.jar ';
	$command .= ' VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true TMP_DIR=/mnt/scratch_mp2/bourque/bourque_group/tmpDir/ ';
	$command .= ' INPUT=' . $laneDirectory . 'alignment/' . $group . '/' . $sampleName . '.sorted.bam ' ;
	$command .= ' OUTPUT=' . $laneDirectory . 'alignment/' . $group . '/' . $sampleName . '.QueryName.bam ' ;
	$command .= ' SORT_ORDER=queryname ';
	
	$retVal{'command'} = $command;
	return ( \%retVal );
}

sub refGtf2matrix {
        my $rH_cfg      = shift;
	my $refGtf  = shift;
	my $readCountDir = shift;
	my $readcountExtension = shift;
	my $outputDir = shift;
	my $outputMatrix  = shift;
	
        my $command = ;
        $command .= 'awk \' BEGIN {FS=\";\"} { split($1,ens,\"\\"\"); split($4,na,\"\\"\") ; print ens[2] \"\t\" na[2]} \'' ;
        $command .= ' ' .$refGtf ;
        $command .= ' | sort -u > ' .$outputDir .'/tmpMatrix.txt ;';
        $command .= ' for i in \` ls ' .$readCountDir .'/*' .$readcountExtension .' \` ;';
        $command .= ' do sort -k1,1 \$i > ' .$outputDir .'/tmpSort.txt ;';
        $command .= ' join -1 1 -2 1' .$outputDir .'/tmpMatrix.txt ' .$outputDir .'/tmpSort.txt > ' .$outputDir .'/tmpMatrix.txt ;';
        $command .= ' done ; rm ' .$outputDir .'/tmpSort.txt ;';
        $command .= ' mv ' .$outputDir .'/tmpMatrix.txt ' .$outputDir .'/' .$outputMatrix .' ;' ;
        $command .= ' rm ' .$outputDir .'/tmpMatrix.txt ';
        
        return $command;
}

1;










































