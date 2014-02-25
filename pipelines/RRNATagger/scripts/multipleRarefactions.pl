#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
itaggerMultipleRarefactions.pl

PURPOSE:
OTU tables rarefaction wrapper. Takes in input multiple OTU tables and
rarefy it to the lowest value of all the OTU tables. Implements the
single_rarefaction.py from Qiime.

INPUT:
--infile <string>          : OTU table in Biom format.
                             Can be multiple infile.
--infileTsv <string>       : OTU table in Qiime (tsv) format.
                             Can be multiple infile.
--threshold <float>        : Fraction of average at which a sample
                             is considered failed. 
--n <int>                  : If table is to be rarefied at <int> reads. 
--step <int>               : Incrementation values in performing rarefactions. 
--permutations <int>       : Number of permutations. Default = 10.
--num_threads <int>        : Num threads.

OUTPUT:
--outdir <string>          : Filtered otu table in Qiime format
                             Must be the same number if outfile
                             as infiles.

NOTES:
Please make sure that the two header lines of the OTU table is 
properly formatted. For instance:

#Full OTU Counts
#OTU ID 01  03  04  05  06  07  08  09  10  11  Consensus lineage	

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $infileTsv, $outdir, $step, $begin, $end, $threshold, $n, $perm, $num_threads);
my $verbose = 0;

GetOptions(
    'infile=s' 		 => \$infile,
	'infileTsv=s'    => \$infileTsv,
	'outdir=s' 		 => \$outdir,
	'threshold=f' 	 => \$threshold,
	'n=i'			 => \$n,
	'permutations=i' => \$perm,
	'num_threads=i'  => \$num_threads,
	'step=i'		 => \$step,
    'verbose' 		 => \$verbose,
    'help' 			 => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile file required\n") unless $infile;
die("--outdir outdir required\n") unless $outdir;
die("--step step required\n") unless $step;
die("--n n required\n") unless $n;
$perm = 10 unless($perm);
$num_threads = 1 unless($num_threads);

my $multipleRarefactions = `which multiple_rarefactions.py`; die if(!defined($multipleRarefactions)); chomp $multipleRarefactions;

$threshold = 0.05 unless($threshold);

## MAIN

# first pass through all infiles to find the lowest abundance value.
my $max = 0;
my $min = 999999999;
my @final_sums;
open(IN, "<".$infileTsv) or die "Can't open infile ".$infileTsv."\n";
my @sums;
my $counter = 0;	

while(<IN>){
	chomp;
	if($_ =~ m/#/){
		next;
	}
	
	my @row = split(/\t/, $_);
	pop(@row); shift(@row);

	# For first row
	if($counter == 0){
		for(my $i=0; $i<@row; $i++){
			$sums[$i] = $row[$i];
		}	
			
		$counter++;
		next;
	}

	for(my $i=0; $i<@row; $i++){
		if(defined($sums[$i])){
			$sums[$i] = $sums[$i] + $row[$i];		
		}else{
			print "Not defined...at indice ".$i." file: ".$infile." line ".$.."\n";
			#die ?
			#$sums[$i] = 0;
			#$sums[$i] = $sums[$i] + $row[$i];		
		}
	}

	$counter++;
}
#Push sums of these samples into the final sums array.
foreach my $sum (@sums){
#	print $sum."\t";
	push(@final_sums, $sum);
}
#print "========================\n";
close(IN);

#For Debug...
foreach(@final_sums){
	#print STDERR "[DEBUG] ".$_."\n";
}

# pick minimal value.
# Find a way to not pick min value if that min value is an outlier.
# For instance, a sample that failed and just have a total of 4 reads
# when the average of sample abundance is say 1500 reads.
@final_sums = sort @final_sums;
my $stat = Statistics::Descriptive::Sparse->new();
foreach my $value (@final_sums){
	$stat->add_data($value);
}
my $avg = $stat->mean();

# Then only consider a minimal value to be good if it is above at least 10% of the average sample abundance.
# This is to reject failed samples.
foreach my $value (@final_sums){
	$min = $value if( ($value < $min) && ($value > $avg * $threshold) );
}
print STDERR "Performing rarefaction at baseline <".$min."> for each tables.\n" if($verbose);		

# Then, normalize all tables upon that value.
# Run single_rarefactions.py 
my $cmd = "parallel_multiple_rarefactions.py";
$cmd .= " -i ".$infile;
$cmd .= " -m 1";
$cmd .= " -x ".$min;
$cmd .= " -s ".$step;
$cmd .= " -n ".$perm;
$cmd .= " -o ".$outdir."/";
$cmd .= " -O ".$num_threads;
print STDERR $cmd."\n";
system($cmd);
$? != 0 ? die "command failed: $!\n" : print STDERR $cmd." done!\n"; 

