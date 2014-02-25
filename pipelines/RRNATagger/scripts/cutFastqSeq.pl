#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
CutFastqSeq.pl

PURPOSE:
To remove a specific number of nucleotides from either
the beginning or end of reads (or both). Useful if we 
have a library that consistently shows bad quality after 
nt 125..., or a library in which the last nucleotides 
are consistently of bad quality.
Also useful to trimm reads before assembly with FLASH.
INPUT:
--infile <string>  : fastq infile
--begin <int>      : Number of nucleotides to remove at the 5 prime region. 
                     Ex: if 3 is entered, the first 3 nucleotides will be removed
--end <int>        : Number of nucleotides to remove at the 3 prime region. 
                     Ex: if 3 is entered, the last 3 nucleotides will be removed
OUTPUT:
--outfile <string> : trimmed fastq file.
NOTES:

BUGS/LIMITATIONS:
TODO Parallelize this script.

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $begin, $end);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'infile=s' 		=> \$infile,
	'outfile=s' 	=> \$outfile,
	'begin=s' 		=> \$begin,
	'end=s' 		=> \$end,
    'verbose' 		=> \$verbose,
    'help' 			=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--fastq fastq file required\n") 			unless $infile;
die("--begin or --end int value required\n") 	unless(defined($begin) or defined($end));
die("--outfile outfile required\n") 			unless $outfile;

## MAIN
open(OUT, ">".$outfile) or die "Can't open file ".$!."\n";

my $in = new Iterator::FastqDb($infile) or die("Unable to open Fastq file, $infile\n");

if(!$begin){
	$begin = 0;
}
if(!$end){
	$end = 0;
}

while( my $curr = $in->next_seq() ){
	my $length = length($curr->seq()) - $end - $begin;
	die "First nucleotides to cut must be higher or equal than/to read length.\n" if($begin >= length($curr->seq()));

	print OUT $curr->header."\n".substr($curr->seq, $begin, $length)."\n+\n".substr($curr->qual, $begin, $length)."\n";
}
close(OUT);

exit;
