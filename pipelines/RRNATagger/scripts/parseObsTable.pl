#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use List::Util qw(sum);

my $usage=<<'ENDHERE';
NAME:
ParseObsTable.pl

PURPOSE:
Script that removes singlets or doublets from an obs table
coming from SeqObs.

INPUT:
Specify one:
--infile_tab <string>     :  SeqObs tab file.
--infile_fasta <string>   :  SeqObs fasta file
--singlet_threshold <int> :  will remove lines with more than <int> singlets 
--doublet_threshold <int> :  will remove lines with more than <int> doublets
--min_value <int>         :
--row_sum  <int>          :  set this argument if rows are to be selected 
                             based on cluster sums. 
OUTPUT:
--outfile_tab <string>    :  Parsed SeqObs table.
--outfile_fasta <string>  :  Parsed SeqObs fasta file.

NOTES:
Lines having at least one column having a value of more 
than 2 (so 3 or more) will be kept.
BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile_tab, $infile_fasta, $singlet_threshold, $doublet_threshold, $outfile_tab, $outfile_fasta, $min_value, $row_sum);
my $verbose = 0;

GetOptions(
    'infile_tab=s' 			=> \$infile_tab,
    'infile_fasta=s' 		=> \$infile_fasta,
    'outfile_tab=s' 		=> \$outfile_tab,
    'outfile_fasta=s' 		=> \$outfile_fasta,
	'min_value=i' 			=> \$min_value,
	'singlet_threshold=i' 	=> \$singlet_threshold,
	'doublet_threshold=i' 	=> \$doublet_threshold,
	'row_sum' 				=> \$row_sum,
    'verbose' 				=> \$verbose,
    'help' 					=> \$help
);
if ($help) { print $usage; exit; }
## VALIDATION

## MAIN
my %hash = ();

open(IN_TAB, $infile_tab) or die "Can't open file ".$infile_tab." ".$!."\n";
open(OUT_TAB, ">".$outfile_tab) or die "Can't open file ".$outfile_tab." ".$!."\n";
open(OUT_FASTA, ">".$outfile_fasta) or die "Can't open file ".$outfile_fasta." ".$!."\n";

if($row_sum){

	my $threshold_value;
	if($min_value){
		$threshold_value = $min_value;
	}else{
		$threshold_value = 1;
	}

	while(<IN_TAB>){
		my $singlets = 0;
		my $doublets = 0;
		my $high_values = 0;
		my $zero = 0;
	
		chomp($_);
	
		if($_ =~ m/\#CLUSTER/){
			print OUT_TAB $_."\n";
			next;
		}
	
		my @row = split("\t", $_);
		my $id = shift(@row);
		my $sum = sum(@row);
		$id = ">".$id;
		#print $id."\n";	
	
		if($sum >= $threshold_value){
			print OUT_TAB $_."\n";
			$hash{$id} = $id;
		}else{

		}	
	}

}else{

	my $threshold_value;
	if($min_value){
		$threshold_value = $min_value;
	}else{
		$threshold_value = 1;
	}
	
	while(<IN_TAB>){
		my $singlets = 0;
		my $doublets = 0;
		my $high_values = 0;
		my $zero = 0;
	
		chomp($_);
	
		if($_ =~ m/\#CLUSTER/){
			print OUT_TAB $_."\n";
			next;
		}
	
		my @row = split("\t", $_);
		my $id = shift(@row);
		$id = ">".$id;
	
		foreach my $el (@row){
			#print $el."\t";
			if($el < $threshold_value){
				$zero++;
			}elsif($el == $threshold_value){
				$singlets++;
			}else{
				$high_values++;
			}
		}
	
		if($high_values >= 1){
			print OUT_TAB $_."\n";
			$hash{$id} = $id;
		#if singlets higher than threshold value, keep row
		}elsif($singlets > $singlet_threshold){
			print OUT_TAB $_."\n";
			$hash{$id} = $id;
		#If not, do not keep it.
		}else{
	
		}
	}
}

#Generate parsed fasta file
my $in = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
while(my $seq = $in->next_seq()){
	if(exists $hash{$seq->header()}){
		print OUT_FASTA $seq->header()."\n".$seq->seq()."\n";
	}
}
exit;
