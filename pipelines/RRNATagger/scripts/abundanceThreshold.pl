#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
AbundanceThreshold.pl

PURPOSE:
From one OTU table, generate new OTU tables having OTUs above certain abundance
threshold only. Abundance threshold calculated will be 10e-6, 10e-5, 10e-4,
10e-3, 10e-2 and 10e-1

Then, will generate a table having different alpha divesity values for each 
abundance threshold values.

INPUT:
--infile <string> : OTU table in Qiime format (qiime.org)
				
OUTPUT:
--outdir <string> : OTU table in Qiime format having OTUs above a certain threshold.
                    Will also contain a spreadsheet with alpha diversity values
                    for each abundance thresholds.

NOTES:


BUGS/LIMITATIONS:
Executing this scripts requires that you have all the required variable, paths, etc. related to Qiime
loaded in your environment. If not it will crash when attempting to compute alpha diversity.
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $outdir);
my $verbose = 0;

GetOptions(
    'infile=s' 	=> \$infile,
	'outdir=s'	=> \$outdir,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die "Provide an --infile arg!\n" unless($infile);
die "Provide an --outdir arg!\n" unless($outdir);

## MAIN
my @thresholds = (10e-6, 10e-5, 10e-4, 10e-3, 10e-2, 10e-1);
my $outdir_data = $outdir."/abundance_thresholds/";
mkdir $outdir_data unless -d $outdir_data;

## Load OTUs abundance into hash
my %hash = ();
my $sum = 0;
my $header = "";
my $prefix = basename($infile);
open(IN, $infile) or die "Can't open file ".$infile."\n";
while(<IN>){
	chomp;
	if($_ =~ m/#/){
		$header .= $_."\n";
		next; 
	}
	my @row = split(/\t/, $_);
	my $id = shift(@row);
	my $lineage = pop(@row);
	my $curr_sum = 0;
	foreach my $el(@row){
		$curr_sum += $el;	
	}
	$hash{$_} = $curr_sum;
	$sum = $sum + $curr_sum;
}

## Generates new OTU tables based on abundance thresholds
## I.e. Just keep OTUs having abundance higher than xx%
my @new_otu_tables;
foreach my $threshold (@thresholds){
	open(OUT, ">".$outdir_data."/".$prefix."_".$threshold) or die "Can't open file ".$outdir_data."/".$prefix."_".$threshold."\n";
	push(@new_otu_tables, $outdir_data."/".$prefix."_".$threshold);
	
	my $new_header = "";

	#$new_header =~ s/OTU ID\t(\S+)/OTU ID\t$1_$threshold/;
	my @header = split(/\t/, $header);
	my $first = shift(@header);
	$new_header .= $first."\t";
	my $last = pop(@header);
	foreach my $el(@header){
		$new_header .= $el."_".$threshold."\t";
	}
	$new_header .= $last;

	print OUT $new_header;

	# For each rows, calculate if abundance
	for my $key (keys %hash){
		my $fraction = ($hash{$key} / $sum);
		#print "VALUE:\t".$hash{$key}."\n";
		#print "SUM:\t".$sum."\n";
		#print "FRACTION:\t".$fraction."\n";
		#print "THRESHOLD:\t".$threshold."\n";
		print OUT $key."\n" if($fraction > $threshold);
	}

	close(OUT);
}

## With these new OTU tables, generate alpha diversity tables.
my $alpha_diversity = "alpha_diversity.py";
#print $outdir_data."/alpha_div\n";
mkdir $outdir_data."/alpha_div" unless -d $outdir_data."/alpha_div";# or die "Could not create dir\n";
my @alpha_div_tables;

foreach(@new_otu_tables){

	## Check if OTU table is empty.
	open(CHECK, $_) or die "Can't open file ".$_."\n";
	my @fh = <CHECK>;
	next if @fh < 3; #Header in OTU tables are on 2 lines.
	close(CHECK);	

	my $prefix = basename($_);
	my $alpha_div_string = $alpha_diversity." ";
	$alpha_div_string .= "-i ".$_." ";
	$alpha_div_string .= "-o ".$outdir_data."/alpha_div/".$prefix." ";
	$alpha_div_string .= "-m chao1,observed_species,shannon,simpson ";
	#print"Computing alpha diversity metrics: ".$alpha_div_string."\n";
	#`$alpha_div_string >/dev/null 2>&1`;
	system($alpha_div_string);
	#$? != 0 ? die "command failed: $!\n" :  print "Command ".$alpha_div_string." successfuly executed\n";
	push(@alpha_div_tables, $outdir_data."/alpha_div/".$prefix);
}

## Finally, generate one single file for all diversity metrics
## Each tables contains one header and one entry (2 lines).
## But header has no # char...
open(OUT, ">".$outdir."/abundance_thresholds.tab") or die "Can't open file ".$outdir."/abundance_thresholds.tab";
print OUT "#Sample\tchao1\tobserved_species\tshannon\tsimpson\n";
foreach(@alpha_div_tables){
	open(ALPHA, $_) or die "Can't open file ".$_."\n";
	my @fh = <ALPHA>;
	my $i=0;
	foreach my $el (@fh){
		print OUT $el  if($i != 0);
		$i++;
	}
	print OUT "\n";
}
close(OUT);
exit;
