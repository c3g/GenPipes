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
use SAMtools;
use File::Basename;

our $rH_cfg;
our $sampleName;
our $rH_laneInfo;
our $readFile;
our $group;

sub matrixMake {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    my $db = shift;    # blast database
    $group       = shift;
    
	$group = ( !defined $group ) ? $sampleName : $group;
    

    # option used if more than one db was specified on the config file.
    # In this case $db should be passed as an argument
    #------------------------------------------------------------------
    $rH_cfg->{'blast.db'} = defined($db) ? $db : $rH_cfg->{'blast.db'};

    my %retVal;
    my $command       = '';
    my $laneDirectory = "read_count/" . $group . "/";

    $command .= ' sh ' . $rH_cfg->{'htseq.tempMatrix'} . ' ' . 'alignment/' . $group . '/' . $group . '.gtf';
    $command .= ' ' . 'assembly/' . $group . '/fasta_split/' . basename($rH_cfg->{'blast.db'}) . '/blast_BestHit.txt ';
    $command .= ' ' . $laneDirectory . 'tmpmatrix.csv ;';
    $command .= ' sh ' . $rH_cfg->{'htseq.fullMatrix'} . ' ' . $laneDirectory . ' ;';
    $command .= ' cp ' . $laneDirectory . 'tmpmatrix.csv DGE/' . $group . '/matrix.csv ;';

    $retVal{'command'} = $command;
    return ( \%retVal );
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


sub sortRead {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $group       = shift;
	
	$group = ( !defined $group ) ? $sampleName : $group;
    my %retVal;
    my $command       = '';
    my $laneDirectory = "alignment/" . $group . "/";

    $command .= ' mkdir -p  read_count/' . $group . ' ;';
    $command .= ' mkdir -p  DGE/' . $group . ';';
    $command .= ' module add '.LoadConfig::getParam($rH_cfg, 'htseq', 'moduleVersion.java') . ' ' . LoadConfig::getParam($rH_cfg, 'htseq', 'moduleVersion.picard') .' ;';
    $command .= ' java -Xmx30g -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'htseq', 'tmpDir').' -jar ${PICARD_HOME}/SortSam.jar ';
    $command .= ' VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true TMP_DIR='.LoadConfig::getParam($rH_cfg, 'htseq', 'tmpDir');
    $command .= ' INPUT=' . $laneDirectory . $sampleName . '.sorted.bam ';
    $command .= ' OUTPUT=read_count/' . $group . '/' . $sampleName . '.QueryName.bam ';
    $command .= ' SORT_ORDER=queryname ';

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
		$strandInfo='no';
	}
	
	my $command ;
	my $latestBam = -M $inputBam;
	my $output1 = -M $outputFile;
	if(!defined($latestBam) || !defined($output1) || $latestBam < $output1) {
		$command .= ' module load ' .LoadConfig::getParam($rH_cfg, 'htseq','moduleVersion.htseq') .' ; ';
		$command .= ' ' .SAMtools::viewFilter($rH_cfg, $inputBam) ;
		$command .= ' | htseq-count - ' .  $inputGtf ;
		$command .= ' -s ' .$strandInfo;
		$command .= ' >' . $outputFile ;
	}
	return $command;
}


sub refGtf2matrix {
        my $rH_cfg      = shift;
	my $refGtf  = shift;
	my $readCountDir = shift;
	my $readcountExtension = shift;
	my $outputDir = shift;
	my $outputMatrix  = shift;
	
        my $command ;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'htseq','moduleVersion.tools') .' &&';
	$command .= ' gtf2tmpMatrix.awk ' .$refGtf;
	$command .= ' ' .$outputDir .'/tmpMatrix.txt &&';
	$command .= ' HEAD=\"Gene\tSymbol\" &&';
	$command .= ' for i in \` ls ' .$readCountDir .'/*' .$readcountExtension .' \` ;';
	$command .= ' do sort -k1,1 \$i > ' .$outputDir .'/tmpSort.txt ;';
	$command .= ' join -1 1 -2 1 ' .$outputDir .'/tmpMatrix.txt ' .$outputDir .'/tmpSort.txt > ' .$outputDir .'/tmpMatrix.2.txt ;';
	$command .= ' mv ' .$outputDir .'/tmpMatrix.2.txt ' .$outputDir .'/tmpMatrix.txt ;';
	$command .= ' na=\$(basename \$i | cut -d\. -f1) ;';
	$command .= ' HEAD=\"\$HEAD\t\$na\" ;';
	$command .= ' done &&';
	$command .= ' echo -e \$HEAD | cat - ' .$outputDir .'/tmpMatrix.txt | tr \' \' \'\t\' > ' .$outputDir .'/' .$outputMatrix .' &&';
	$command .= ' rm ' .$outputDir .'/tmpSort.txt ' .$outputDir .'/tmpMatrix.txt ';
        
        return $command;
}

1;

