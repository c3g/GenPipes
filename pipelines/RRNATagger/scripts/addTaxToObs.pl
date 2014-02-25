#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $usage=<<'ENDHERE';
NAME:
itaggerAddTaxToSeqobs.pl

PURPOSE:
This script takes a SeqObs abundance table and it's corresponding 
RDP taxonomy (generated with the RDP classifier tool).
It appends the RDP lineages to each corresponding SeqObs clusters 
and format the output as a OTU table compatible with the Qiime suite tools.

INPUT:
-seqobs <tab_file>         :  SeqObs output tablablature file (cluster abundance/indexes)
-rdp <tab_file>            :  RDP output
-tax_level <string>        :  Taxonomic level to apply selection cutoff, either "species", "genus"
                              , "family", "order", "class", "phylum" or "best"without the quotes (""). 
                              If this argument is not set, will return the deepest node in the lineage 
                              that meets the threshold value.			
-cutoff <int>              :  Selection cutoff, value between 0.0 and 1. Default is 0.8.

OUTPUT:
--outfile <tab_file>        :  OTU table compatible with the Qiime tools.
--outfile_failed <tab_file> :  OTU table compatible with the Qiime tools. 
                               Contains values that didn't meet the cutoff requirement.
NOTES:
If --tax_level arg is set to "best", lineages of different lengths
will be returned. Will return the lineage all the way to the root 
to the deepest node that pass the cutoff value.

BUGS/LIMITATIONS:
	
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $seqobs, $rdp, $outfile, $outfile_failed, $tax_level, $cutoff);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'seqobs=s' => \$seqobs,
	'rdp=s' => \$rdp,
	'tax_level=s' => \$tax_level,
	'cutoff=f' => \$cutoff,
	'outfile=s' => \$outfile,	
	'outfile_failed=s' => \$outfile_failed,	
    'verbose' => \$verbose,
    'help' => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--seqobs seqobs table (qiime format) file required.\n") unless $seqobs;
die("--rdp rdp table required.\n") unless $rdp;
die("--outfile outfile required.\n") unless $outfile;
$cutoff = 0.80 unless $cutoff;
#$tax_level = 0 unless $tax_level;
die("--cutoff cutoff value must be between 0 and 1.0 inclusively.\n") unless($cutoff >= 0.00 or $cutoff <= 1.00);
die("--tax_level tax_level is not an acceptable value. See --help for info.\n") unless($tax_level ne "species" or $tax_level ne "genus" or $tax_level ne "family" or $tax_level ne "order" or $tax_level ne "class" or $tax_level ne "phylum" or $tax_level ne "kingdom" or $tax_level ne "best");
$cutoff = 0.00 if($cutoff == 0 or $cutoff == 0.0);
#$cutoff = int($cutoff);

## MAIN
print STDERR "Appending RDP taxonomy to SeqObs abundance table...\n" if $verbose == 1;

open(OUT, ">".$outfile) or die "Can't open file ".$outfile."\n";
open(OUT_FAILED, ">".$outfile_failed) or die "Can't open file ".$outfile_failed."\n";
open(SEQOBS, $seqobs) or die "Can't open file ".$seqobs."\n";
open(RDP, $rdp) or die "Can't open file ".$rdp."\n";

print OUT "#Full OTU Counts\n";
print OUT_FAILED "#Full OTU Counts\n";

my %hash = ();
my %hash_prob = ();
my %hash_prob_best = ();

#Put relevant rdp entries in a hash table
while(<RDP>){
	chomp($_);
	my @row = split(/\t/, $_);
	#print $_."\n" if($verbose);
	foreach my $el (@row){
		#print $el."\n" if($verbose);
	}
		
	#check if entered tax_level is concordant with the actual number of tax levels in the RDP file
	if($tax_level ne "best"){
		die "This RDP file contains too many taxonomic levels, this script supports 7 taxonomic levels.\n" if(@row > 26 );
		die "This RDP file's lineage does not go include the species level.\n" if($tax_level eq "species" && @row < 26);
		die "This RDP file's lineage does not go include the genus level.\n" if($tax_level eq "genus" && @row < 23);
		die "This RDP file's lineage does not go include the family level.\n" if($tax_level eq "family" && @row < 20);
		die "This RDP file's lineage does not go include the order level.\n" if($tax_level eq "order" && @row < 17);
		die "This RDP file's lineage does not go include the class level.\n" if($tax_level eq "class" && @row < 14);
		die "This RDP file's lineage does not go include the phylum level.\n" if($tax_level eq "phylum" && @row < 11);
		die "This RDP file's lineage does not go include the kingdom level.\n" if($tax_level eq "kingdom" && @row < 8);
	}
		#Output of RDP classifier guideline...
		#$row[5]: Kingdom
		#$row[7]: k_value
		#$row[8]: Phylum
		#$row[10]: p_value
		#$row[11]: Class
		#$row[13]: c_value
		#$row[14]: Order
		#$row[16]: o_value
		#$row[17]: Family
		#$row[19]: f_value
		#$row[20]: Genus
		#$row[22]: g_value
		#$row[23]: Species
		#$row[25]: s_value
	
	my $id = $row[0];
	$id =~ s/@//;
	$id =~ s/>//;
	if(substr($id, -1) eq "#"){
		$id = substr($id, 0, -1);
	}
	$id = $1 if($id =~ m/(\d+)/);
	
	print STDERR "[DEBUG]: ".$id."\n" if($verbose);

	if(@row == 26){
		$hash{$id} = $row[5].";".$row[8].";".$row[11].";".$row[14].";".$row[17].";".$row[20].";".$row[23];
		$hash_prob_best{$id} = $row[7]."\t".$row[10]."\t".$row[13]."\t".$row[16]."\t".$row[19]."\t".$row[22]."\t".$row[25];
		#print STDERR $hash_prob_best{$id}."\n" if($verbose);
	}
	if(@row == 23){
		$hash{$id} = $row[5].";".$row[8].";".$row[11].";".$row[14].";".$row[17].";".$row[20];
		$hash_prob_best{$id} = $row[7]."\t".$row[10]."\t".$row[13]."\t".$row[16]."\t".$row[19]."\t".$row[22]."\t"."0";
		#print $id."\t".$row[5].";".$row[8].";".$row[11].";".$row[14].";".$row[17].";".$row[20]."\n";
	}
	if(@row == 20){
		$hash{$id} = $row[5].";".$row[8].";".$row[11].";".$row[14].";".$row[17];
		$hash_prob_best{$id} = $row[7]."\t".$row[10]."\t".$row[13]."\t".$row[16]."\t".$row[19]."\t"."0"."\t"."0";
	}
	if(@row == 17){
		$hash{$id} = $row[5].";".$row[8].";".$row[11].";".$row[14];
		$hash_prob_best{$id} = $row[7]."\t".$row[10]."\t".$row[13]."\t".$row[16]."\t"."0"."\t"."0"."\t"."0";
	}
	if(@row == 14){
		$hash{$id} = $row[5].";".$row[8].";".$row[11];
		$hash_prob_best{$id} = $row[7]."\t".$row[10]."\t".$row[13]."\t"."0"."\t"."0"."\t"."0"."\t"."0";
	}
	if(@row == 11){
		$hash{$id} = $row[5].";".$row[8];
		$hash_prob_best{$id} = $row[7]."\t".$row[10]."\t"."0"."\t"."0"."\t"."0"."\t"."0"."\t"."0";
	}
	if(@row == 8){
		$hash{$id} = $row[5];
		$hash_prob_best{$id} = $row[7]."\t"."0"."\t"."0"."\t"."0"."\t"."0"."\t"."0"."\t"."0";
	}
	
	#print "ROW: ".@row."\n" if($verbose);

	if($tax_level ne "best"){
		$hash_prob{$id} = ($row[7]) if($tax_level eq "kingdom");	
		$hash_prob{$id} = ($row[10]) if($tax_level eq "phylum");	
		$hash_prob{$id} = ($row[13]) if($tax_level eq "class");	
		$hash_prob{$id} = ($row[16]) if($tax_level eq "order");	
		$hash_prob{$id} = ($row[19]) if($tax_level eq "family");	
		$hash_prob{$id} = ($row[22]) if($tax_level eq "genus");	
		$hash_prob{$id} = ($row[25]) if($tax_level eq "species");
	}#else{ #Keep the deepest level/node that meets threshold

	#}
}
close(RDP);

#Loop through seqobs table and add rdp taxonomy accordingly
my $i=0; #for the header.
while(<SEQOBS>){
	chomp($_);
	if($i == 0){ #Print Qiime-style OTU table header.
		$_ =~ s/CLUSTER/OTU ID/;
		print OUT $_."\ttaxonomy\n"; 
		print OUT_FAILED $_."\ttaxonomy\n"; 
	}else{
		my @row = split(/\t/, $_);
		my @tax;
		
		my $id = $row[0];
		$id =~ s/@//;
		$id =~ s/>//;
		if(substr($id, -1) eq "#"){
			$id = substr($id, 0, -1);
		}
		print $id."\n" if($verbose);
		if(exists ($hash{$id}) ){
			@tax = split(/;/, $hash{$id});
		}else{
			print STDERR "Undefined taxonomy for : ".$id."==> may be cause by reads/clusters shorter than 50 nt.\n";
		}
	
	
		if(defined $hash_prob{$id} or defined $hash_prob_best{$id}){
			if($tax_level ne "best"){
				my $tax;	
				$tax = $tax[0].";".$tax[1].";".$tax[2].";".$tax[3].";".$tax[4].";".$tax[5].";".$tax[6] if($tax_level eq "species");
				$tax = $tax[0].";".$tax[1].";".$tax[2].";".$tax[3].";".$tax[4].";".$tax[5] if($tax_level eq "genus");
				$tax = $tax[0].";".$tax[1].";".$tax[2].";".$tax[3].";".$tax[4] if($tax_level eq "family");
				$tax = $tax[0].";".$tax[1].";".$tax[2].";".$tax[3] if($tax_level eq "order");
				$tax = $tax[0].";".$tax[1].";".$tax[2] if($tax_level eq "class");
				$tax = $tax[0].";".$tax[1] if($tax_level eq "phylum");
				$tax = $tax[0] if($tax_level eq "kingdom");

				#print $hash_prob{$row[0]}."\t".$cutoff."\n";
				($hash_prob{$id} >= $cutoff)? print OUT $_."\t".$tax."\n" : print OUT_FAILED $_."\t".$hash{$id}."\n";	
				#add step to check if hash is empty?... could be possible?
			}else{
				my @probs = split(/\t/, $hash_prob_best{$id});

				my $indice=0;
				my $tax="";
				while (@probs and $indice < @tax ){
					#print $probs[$indice]."\t";
					if($probs[$indice] >= $cutoff){
						$tax .= $tax[$indice].";"
					}else{
						last;
					}
					$indice++;
				}
				#print "\n";
				print OUT $_."\t".$tax."\n" if($indice > 0);
				print OUT_FAILED $_."\t".$hash{$id}."\n" if($indice == 0);
			}
		}else{
			#It can happen that if read length is too short (eg, < 50 nt), it wont be classified by RDP.
		}
	}
	$i=1;
}
close(SEQOBS);

exit;

