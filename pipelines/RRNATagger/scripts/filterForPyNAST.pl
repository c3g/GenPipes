#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
filterForPyNAST.pl

PURPOSE:

INPUT:
--infile_otu_table <string>     : OTU table in Qiime format (qiime.org)
--infile_fasta <string>         : Fasta file
				
OUTPUT:
STDOUT

NOTES:


BUGS/LIMITATIONS:

 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile_otu_table, $infile_fasta);
my $verbose = 0;

GetOptions(
    'infile_otu_table=s'  => \$infile_otu_table,
	'infile_fasta=s'      => \$infile_fasta,
    'verbose' 	          => \$verbose,
    'help' 		          => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "--infile_otu_table <string> required\n" unless($infile_otu_table);
die "--infile_fasta <string> required\n" unless($infile_fasta);

## MAIN
my %hash_otu_table;
open(OTU_IN, "<".$infile_otu_table) or die "Can't open OTU table ".$infile_otu_table."\n";
while(<OTU_IN>){
	chomp;
	next if($_ =~ m/#/);
	my @row = split(/\t/, $_);
	$hash_otu_table{$row[0]} = $_;	
}
close(OTU_IN);

my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $header = $curr->header();
	$header =~ s/>//;
	if(exists $hash_otu_table{$header}){
		print STDOUT ">".$header."\n".$curr->seq."\n";
	}
}
exit;
