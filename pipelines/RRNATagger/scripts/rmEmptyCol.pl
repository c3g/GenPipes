#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
RmEmptyCol.pl

PURPOSE:
Removes empty columns in a otu table. This step is needed to generate 
graphs with Qiime tools. If you have a otu table, but with a columns 
having only zeros, taxa summary will not be plotted.

INPUT:
--infile <fastq_infile>    : otu table in Qiime format

OUTPUT:
--outfile <fasta_outfile>  : filtered otu table in Qiime format

NOTES:
Please make sure that the two header lines of the OTU table is 
properly formatted. For instance:

#Full OTU Counts
#OTU ID 01  03  04  05  06  07  08  09  10  11  Consensus lineage	

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $begin, $end);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'infile=s' 	=> \$infile,
	'outfile=s' => \$outfile,
	'begin=s' 	=> \$begin,
	'end=s' 	=> \$end,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

#VALIDATE
die("--infile file required\n") unless $infile;
die("--outfile outfile required\n") unless $outfile;

#MAIN=====================================================================================================================================================================
open(IN, $infile) or die "Can't open file ".$!."\n";
open(OUT, ">".$outfile) or die "Can't open file ".$!."\n";

my %hash = ();

#Loop a first time through the otu table.
while(<IN>){
	chomp($_);
	if($_ =~ m/\#/){
		next;
	}
	my @row = split(/\t/, $_);
	shift(@row);
	pop(@row);
	my $i = 0;
	foreach my $el (@row){ #Summarize the columns.
		#$hash{$i} = 0 if(!defined $hash{$i});
		$hash{$i} = 0 if(!exists $hash{$i});
		$hash{$i} = $hash{$i} + $el;
		$i++;
	}	
}
close(IN);

open(IN, $infile) or die "Can't open file ".$!."\n";
while(<IN>){
	chomp($_);
	if($. == 1){
		print OUT $_."\n";
		next;
	}
	
	my @row = split(/\t/, $_);
	print OUT shift(@row)."\t";
	my $lineage = pop(@row);
	my $i = 0;
	foreach my $el(@row){
	print OUT $el."\t" if($hash{$i} != 0);
		$i++;
	}	
	print OUT $lineage."\n";
}
close(IN);

exit;
