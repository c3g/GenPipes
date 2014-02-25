#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
FilterObsTable.pl

PURPOSE:
Keep OTU that have an abundance of at least <threshold> at least <frequency> times.
For instance, keep OTUs that have at least 9 reads on at least 2 samples.

INPUT:
--obs_tab <string>       : Obs tabular file (output from SeqOBs)
--obs_fasta <string>     : Fasta file (output from SeqObs)
OR
--otu_table <string>     : OTU table.
--threshold <int>        : abundance threshold (i.e. 9).
--frequency <int>        : frequency threhsold. (i.e. 1).

OUTPUT:
--obs_tab_out <string>   : Obs tabular file with no singlets
--obs_fasta_out <string> : Fasta file with no singlets sequences.
OR
--otu_table_out <string> : Filtered OTU table.

NOTES
For OTU table filtering. OTUs must have an abundance higher than --threshold 
in at least --frequency samples.

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $obs_tab, $obs_fasta, $obs_tab_out, $obs_fasta_out, $otu_table, $otu_table_out, $threshold, $frequency);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'obs_tab=s' 		=> \$obs_tab,
	'obs_fasta=s' 		=> \$obs_fasta,
	'obs_tab_out=s' 	=> \$obs_tab_out,
	'otu_table=s'		=> \$otu_table,
	'threshold=i'		=> \$threshold,
	'frequency=i'		=> \$frequency,
	'otu_table_out=s'	=> \$otu_table_out,
	'obs_fasta_out=s' 	=> \$obs_fasta_out,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
if($otu_table){
	die("--frequency missing\n") 			unless $frequency;
	die("--threshold missing\n") 			unless $threshold;
	die("--otu_table_out file required\n") 	unless $otu_table_out;
}else{
	die("--obs_tab file required\n")	 	unless $obs_tab;
	die("--obs_fasta file required\n") 		unless $obs_fasta;
	die("--obs_tab_out file required\n") 	unless $obs_tab_out;
	die("--obs_fasta_out file required\n") 	unless $obs_fasta_out;
}

## MAIN
if($otu_table){
	open(IN, "<".$otu_table) or die "Can't open file ".$otu_table."\n";
	open(OUT, ">".$otu_table_out) or die "Can't open file ".$otu_table_out."\n";
	while(<IN>){
		chomp;
		if($_ =~ m/#/){
			print OUT $_."\n";
			next;
		}
		my @row = split(/\t/, $_);
		shift(@row);pop(@row);
		my $count = 0;
		foreach my $el(@row){
			if($el >= $threshold){
				$count++;
			}
		}	
		
		if($count >= $frequency){
			print OUT $_."\n";
		}
	}
	close(IN);
	close(OUT);
	
}else{
	open(IN_TAB, $obs_tab) or die "Can't open file ".$obs_tab."\n";
	open(OUT_TAB, ">".$obs_tab_out) or die "Can't open file ".$obs_tab_out."\n";
	open(OUT_FASTA, ">".$obs_fasta_out) or die "Can't open file ".$obs_fasta_out."\n";
	
	my %hash = ();
	while(<IN_TAB>){
		chomp($_);
		if($. == 1){ #skip header
			print $_."\n";
			next;	
		}
		my @row = split(/\t/, $_);
		my $id = shift(@row); #remove cluster identifiers
		my $sum = eval join '+', @row;
		if($sum <= 1){ #remove singlets
			print OUT_TAB $_. "\n";
			$hash{$id} = $sum; 
		}
	}
	close(IN_TAB);
	close(OUT_TAB);
	
	
	my $in = new Iterator::FastaDb($obs_fasta) or die("Unable to open Fastq file, $obs_fasta\n");
	while( my $curr = $in->next_seq() ){
		if(exists $hash{$curr->header()}){
			print OUT_FASTA $curr->header()."\n".$curr->seq()."\n";
		}
	}
	close(OUT_FASTA);
}
exit;
