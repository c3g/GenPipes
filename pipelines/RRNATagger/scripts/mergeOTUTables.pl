#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp;
use Iterator::FastaDb;
use Env qw/TMPDIR/;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
MergeOTUTables.pl

PURPOSE:
Merge different OTU tables together. For each OTU table, you
have to submit its equivalent fasta cluster representative file.
After the merge, original cluster/OTU IDs will be lost (A new
unique ID will be given to make sure that there is no conflict - 
for instance if two OTU in two different OTU tables have the 
same ID number.

INPUT:
--infile_otu_table <string>  : OTU table in Qiime format.
                               Can be multiple infile.
--infile_fasta <string>      : Fraction of average at which a sample
                               is considered failed. 

OUTPUT:
--outfile_otu_table <string> : Merged OTU table.
--outfile_fasta     <string> : Merge fasta file.

NOTES:
IMPORTANT: This script assumes that you have the exact same clusters in 
both otu_table and fasta files.

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, @infile_otu_table, @infile_fasta, $outfile_otu_table, $outfile_fasta);
my $verbose = 0;

GetOptions(
    'infile_otu_table=s' 	=> \@infile_otu_table,
	'infile_fasta=s'		=> \@infile_fasta,
	'outfile_otu_table=s' 	=> \$outfile_otu_table,
	'outfile_fasta=s'		=> \$outfile_fasta,
    'verbose' 				=> \$verbose,
    'help' 					=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile_otu_table file required\n")  unless @infile_otu_table;
die("--infile_fasta file required\n")      unless @infile_fasta;
die("--outfile_otu_table file required\n") unless $outfile_otu_table;
die("--outfile_fasta file required\n")     unless $outfile_fasta;
die("must be equal number of --infile_otu_table and --infile_fasta args...\n") if(@infile_otu_table != @infile_fasta);

## MAIN

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerMergeOTUTablesXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);
print STDERR $tmpdir."\n" if($verbose);

# Open outfiles
open(FASTA, ">".$outfile_fasta) or die "Can't open file ".$outfile_fasta."\n";

my $i=0;
my $temp_file_counter = 0;
my @otu_tables;

my %hash;

# Loop through each infile_fasta and infile_otu_table and assign new $i number as OTU ids.
foreach my $infile_otu_table (@infile_otu_table){
	my $infile_fasta = shift(@infile_fasta);

	# Then deal with the otu table.
	my $curr_otu_table = $tmpdir."/otu_table_".$temp_file_counter; 
	push(@otu_tables, $curr_otu_table);
	open(CURR_OTU, ">".$curr_otu_table) or die "Can't open file ".$curr_otu_table."\n";
	$temp_file_counter++;

	open(IN, "<".$infile_otu_table) or die "Can't open file ".$infile_otu_table."\n";
	while(<IN>){
		chomp;
		if($_ =~ m/#/){
			print CURR_OTU $_."\n";
		}else{
			my @row = split(/\t/, $_);					
			my $old_id = shift(@row);
			unshift(@row, $i);
			print CURR_OTU join("\t", @row)."\n";
			$hash{$old_id} = $i;
			$i++;
		}
	}
	close(IN);
	close(CURR_OTU);
	
	my $ref_fasta_db = Iterator::FastaDb->new($infile_fasta) or die("Unable to open Fasta file, $infile_fasta\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		my $header = $curr->header();
		my $old_id;
		if($header =~ m/^>(\d+)/){
			$old_id = $1;
		}else{
			die "Can't find regex...\n";
		}

		if(exists $hash{$old_id}){
			#print STDERR $old_id."\t".$hash{$old_id}."\n";
			print FASTA ">".$hash{$old_id}."\n".$curr->seq."\n";
			delete $hash{$old_id};
		}
	}
}
close(FASTA);

# Then run Qiime's merge otu tables scripts.
my $first_otu_table = shift(@otu_tables);
my $cmd = "merge_otu_tables.py -i ".$first_otu_table;
foreach (@otu_tables){
	$cmd .= ",".$_;
}
$cmd .= " -o ".$outfile_otu_table;
system($cmd);

## REMOVE TEMP FILES
sub END{
	#system("rm ".$tmpdir." -rf");
}

exit;
