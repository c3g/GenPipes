#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastqDb;
use Iterator::FastaDb;
use List::Util qw(sum);

my $usage=<<'ENDHERE';
NAME:
DummyObsTable.pl

PURPOSE:
Generates a dummy Seqobs format table.

INPUT:
--fasta <string>   : fasta file
OR
--fastq <string>   : fastq file

OUTPUT:
--outfile <string> : Seqobs format table

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $fastq, $fasta, $outfile, $complete_header);
my $verbose = 0;

## SCRIPTS
GetOptions(
	'fastq=s' 			=> \$fastq,
	'fasta=s' 			=> \$fasta,
	'outfile=s' 		=> \$outfile,
	'complete_header' 	=> \$complete_header,
    'verbose'			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("Enter either --fasta OR --fastq\n") if($fasta && $fastq);
die("--outfile required\n") unless $outfile;

## MAIN
open(OUT, ">".$outfile) or die "Can't find file ".$outfile."\n";
my @array=();

if($fasta){
	my $fasta_db = Iterator::FastaDb->new($fasta) or die("Unable to open Fasta file, $fasta\n");
	while( my $seq = $fasta_db->next_seq() ) {
		my $id = $seq->header();
		$id =~ s/>//;
		push(@array, $id);
	}
}elsif($fastq){
	my $fastq_db = Iterator::FastqDb->new($fastq) or die("Unable to open Fasta file, $fastq\n");
	while( my $seq = $fastq_db->next_seq() ) {
		my $id;
		if(defined($seq->barcode())){
			$id = "@".$seq->base()."#".$seq->barcode();
		}else{
			$id = "@".$seq->base();
		}
		#print $id."\n";
		
		if($complete_header){
			$id = $seq->header();
			#$id =~ s/@//; #Keep the @ because rdp_classifier 2.5 also keeps it.
			if($id =~ m/(\S+)/){
				$id = $1;
				#print $id."\n";
			}else{
				print STDERR "Did not matched regexp\n";
			}
		}

		push(@array, $id);
	}
}

print OUT "#CLUSTER\tcount\n";
my $i=0;
while($i < @array){
	print OUT $array[$i]."\t1\n";
	$i++;
}
