#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
join_paired_fastqs.pl

PURPOSE:

INPUT:
--reads_1 <fastq_file_containing_1st_reads>
--reads_2 <fastq_file_containing_2nd_reads>
		
OUTPUT:
 --outfile <string> : paired fastq file
NOTES:

BUGS/LIMITATIONS:
This script assumes that both fastq files are in the EXACT SAME order! 

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $end1, $end2, $outfile);
my $verbose = 0;

GetOptions(
    'reads_1=s' => \$end1,
    'reads_2=s' => \$end2,
    'outfile=s' => \$outfile,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## MAIN
open(OUT, ">".$outfile) or die "Can't open file ".$outfile."\n";

my $in_1 = new Iterator::FastqDb($end1,{trimN=>0, trim3=>0});
my $in_2 = new Iterator::FastqDb($end2,{trimN=>0, trim3=>0});
while( my $curr_1 = $in_1->next_seq() ){
	my $curr_2 = $in_2->next_seq();
	print OUT "@".$curr_1->base()."#".$curr_1->barcode()."/".$curr_1->pair()."\n".$curr_1->seq()."\n+\n".$curr_1->qual()."\n";
	print OUT "@".$curr_2->base()."#".$curr_2->barcode()."/".$curr_2->pair()."\n".$curr_2->seq()."\n+\n".$curr_2->qual()."\n";
	die "Problem in file seq order...\n" if($curr_1->base() ne $curr_2->base());
}
close(OUT);

