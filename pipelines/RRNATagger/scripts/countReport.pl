#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;

my $usage=<<'ENDHERE';
NAME:
PURPOSE:

INPUT:
--file
--name 
--analysisType <string>    : Name for that report (i.e. assembled reads, reads1 or reads2).
--barcodesDist <tab_file>  : Barcodes distribution
--OTUtable <string>        : Classified OTUs/clusters. (OTU table).
--obsTable <int>           : Obs. table.

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:
	
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, @analysisType, @barcodesDist, @OTUtable, @obsTable, @name, @file);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'analysisType=s' 	=> \@analysisType,
	'barcodesDist=s' 	=> \@barcodesDist,
	'OTUtable=s' 		=> \@OTUtable,
	'obsTable=s'		=> \@obsTable,
    'name=s' 			=> \@name,
	'file=s' 			=> \@file,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE

## SUBS
# Return commat formatted numbers
# $_[0] - a number
# @returns comma-formatted number
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

# Count number of lines in a file (it is assumed it is a fastq file) and divide by 4. Return the value.
# @input: one fastq file
# @output: a number
sub countFastqReads{
	my $file = $_[0];
	my $count = 0;

	#if((!-s $file) and (!-e $file)){
		if($file =~ m/\.gzip|\.gz/) {
			$count = `gunzip -c $file | cut -c -2 | wc -l`;
		}else{
			$count = `cut -c -2 $file | wc -l`;
		}
		chomp($count);
		$count =~ m/(\d+)/;
		$count = $1 / 4;
		$count =~ s/\.\d{2}//;#temporary patch to remove decimal for potential disrupted fastqs.
		$count = commify($count);
	#}else{
	#	die "$file does not exists...\n"
	#}   
	$count = commify($count);
	
	return $count;
}

## MAIN

# Write sequence/reads reports on seqobs output
# In two cluster files 1-clustered 2-classified clusters 3- barcode dist
# @returns: null
foreach my $name (@name){
	my $file = shift(@file);

	print STDOUT $name."\t".countFastqReads($file)."\n";
	
}		


# Write sequence/reads reports on seqobs output
# In two cluster files 1-clustered 2-classified clusters 3- barcode dist
# @returns: null
foreach my $analysisType (@analysisType){
	my $barcodesDist = shift(@barcodesDist);
	my $OTUtable = shift(@OTUtable);
	my $obsTable = shift(@obsTable);
	
	print STDERR "Barcodes distribution\n" if($verbose);
		
	$analysisType =~ s/_//;
	print STDOUT "Cluster counts for $analysisType\n";
	print STDOUT "===BARCODES DISTRIBUTION===\n";
	open(IN, $barcodesDist) or die "Can't open file ".$barcodesDist."\n";
	while(<IN>){
		chomp($_);
		my @row = split(/\t/, $_);
		print STDOUT $row[0]."\t".$row[1]."\t".commify($row[2])."\t".$row[3]."\n";
	}
	close(IN);

	open(IN, $obsTable) or die "Can't open file ".$obsTable."\n";
	my $line_number = 0;
	my $sum=0;
	while(<IN>){
		$line_number++;
		chomp($_);
		next if ($_ =~ m/\#/);
		my @row = split(/\t/, $_);
		shift(@row);
		foreach my $number (@row){
			$sum += $number;
		}
	}
	close(IN);
	print STDOUT "Clustered:\t".commify($sum)." sequences give a total of ".commify($line_number)." clusters\n";
	
	open(IN, $OTUtable) or die "Can't open file ".$OTUtable."\n";
	$line_number = 0;
	$sum=0;
	while(<IN>){
		chomp($_);
		next if ($_ =~ m/\#/);
		$line_number++;
		my @row = split(/\t/, $_);
		shift(@row);
		pop(@row);
		foreach my $number (@row){
			$sum += $number;
		}
	}	
	close(IN);
	print STDOUT "Classifed clusters:\t".commify($sum)." sequences packed in ".commify($line_number)." cluster\n";
}
exit;
