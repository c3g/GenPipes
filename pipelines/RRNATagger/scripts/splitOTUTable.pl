#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
itaggerSplitOTUTable.pl

PURPOSE:

INPUT:
--infile <string>     : OTU table in Qiime format (qiime.org)
--select <string>     : Either 'bactArch' or 'fungi'
				
OUTPUT:
--matched <string>    : OTU table having only Bacteria AND Archeal organisms.
--unmatched <string>  : OTU table having non-Bacteria and Non-Archeal organisms.

NOTES:


BUGS/LIMITATIONS:

 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $matched, $unmatched, $select);
my $verbose = 0;

GetOptions(
    'infile=s' 	   => \$infile,
	'matched=s'    => \$matched,
	'unmatched=s'  => \$unmatched,
	'select=s'     => \$select,
    'verbose' 	   => \$verbose,
    'help' 		   => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--infile <string> required\n" unless($infile);
die "--matched <string> required\n" unless($matched);
die "--unmatched <string> required\n" unless($unmatched);
die "--select must be either 'bactArch' or 'fungi'" if($select ne "bactArch" && $select ne "fungi");

## MAIN

# TODO Put that in a separate script.
# This subroutine takes a otu_table in input and gives a otu table minus every lineage that does not start with k__Bacteria or k__Archaea in output.
# @input : otu_table.tab
# @output : otu_table.tab (only containing bacteria/archaea organisms).
open(IN, "<".$infile) or die "Can't open input file ".$infile."\n";
open(OUT, ">".$matched) or die "Can't open input file ".$matched."\n";
open(OUT_F, ">".$unmatched) or die "Can't open input file ".$unmatched."\n";

my $regex; 
$regex = "k__bacteria|k__archaea" if($select eq "bactArch");
$regex = "k__fungi" if($select eq "fungi");

while(<IN>){
	chomp;
	if($_ =~ m/#/){
		print OUT $_."\n";
		print OUT_F $_."\n";
		next;
	}
	if($_ =~ m/($regex)/i ){
		print OUT $_."\n";
	}else{
		print OUT_F $_."\n";
	}
}
close(IN);
close(OUT);
close(OUT_F);
