#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
get_subset.pl
PURPOSE:

INPUT:
--infile <string>    : Sequence file (fastq). If interleaved paired end reads
OR
--infile_1 <string>   : Sequence file (fastq). If reads 1 and reads 2 are in separate files
--infile_2 <string>   : Sequence file (fasta). If reads 1 and reads 2 are in separate files

--n <int>            : Number of sequences to get from subset.
--paired             : Set this arg if fastq is paired-end ordered
                       For instance: @blabla1/1
                                   ACGTACGAGCGT
                                   +
                                   HHHH%!@%@HHJ
                                   blabla1/2
                                   ACGTACGAGCGT
                                   +
                                   HHHH%!@%@HHJ

OUTPUT:
--outfile <string>   : Sequence file (fastq). If interleaved paired end reads. 
OR
--outfile_1 <string> : Sequence file (fastq). If reads 1 and reads 2 are in separate files
--outfile_2 <string> : Sequence file (fastq). If reads 1 and reads 2 are in separate files

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $n, $outfile, $paired, $infile_1, $infile_2, $outfile_1, $outfile_2);
my $verbose = 0;

GetOptions(
    'infile=s' 	  => \$infile,
	'infile_1=s'  => \$infile_1,
	'infile_2=s'  => \$infile_2,
    'outfile=s'   => \$outfile,
    'outfile_1=s' => \$outfile_1,
    'outfile_2=s' => \$outfile_2,
	'paired' 	  => \$paired,
	'n=i' 		  => \$n,
    'verbose' 	  => \$verbose,
    'help' 		  => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
if($outfile && $infile){
	open(OUT, ">".$outfile) or die "Can't open file ".$outfile."\n";
	die("File $infile does not exist or is empty...\n") if((!-e $infile) and (!-s $infile));	

}elsif( ($outfile_1 && $outfile_2) && ($infile_1 && $infile_2) ){
	open(OUT_1, ">".$outfile_1) or die "Can't open file ".$outfile_1."\n";
	open(OUT_2, ">".$outfile_2) or die "Can't open file ".$outfile_2."\n";
	die("File $infile_1 does not exist or is empty...\n") if((!-e $infile_1) and (!-s $infile_1));	
	die("File $infile_2 does not exist or is empty...\n") if((!-e $infile_2) and (!-s $infile_2));	

}else{
	die("Please specify either --outfile OR --outfile_1 and --outfile_2\n");
}

die("Please specify --paired with --infile and --outfile") if($paired and $infile_1);
die("--n arg missing...\n") unless($n);

if($paired && $infile){
	## GET NUMBER OF ENTRIES IN THE FILE
	my $range = count_fastq_reads($infile);
	$range = $range / 2; #Because paired sequence file.
	print "Number of entries: ".$range."\n";
	
	## GENERATES N RANDOM NUMBERS BETWEEN 0-NUMBER OF ENTRIES
	my %hash;
	for(my $i=0; $i<$n; $i++){
		my $flag = 0;	
		while($flag == 0){
			my $random_number = int(rand($range)) * 2;
			if(exists $hash{$random_number}){
				$flag = 0;
			}else{
				$hash{$random_number} = $random_number;
				$flag = 1;
			}
		}
	}
	
	my $counter = 0;
	foreach(sort{$a <=>$b} keys %hash){
		print $_."\n";
		$counter++;
	}
	print "Numbers generated: ".$counter."\n";
	
	
	$counter = 0;
	my $ref_fastq_db = Iterator::FastqDb->new($infile, {trimN=>0, trim3=>0}) or die("Unable to open Fastq file, $infile\n");
	my $flag = 0;
	while( my $curr = $ref_fastq_db->next_seq() ) {
		if($flag == 1){
			print OUT $curr->output;
			$flag = 0;
		}

		if(exists $hash{$counter}){
			print OUT $curr->output;
			$flag = 1;
		}
		$counter++;
	}
}elsif( $infile || $infile_1 ){
	## GET NUMBER OF ENTRIES IN THE FILE
	my $curr_infile;
	$curr_infile = $infile if($infile);
	$curr_infile = $infile_1 if($infile_1);
	my $range = count_fastq_reads($curr_infile);
	print STDERR "Number of entries: ".$range."\n" if($verbose);
	$n =  $range if($range <= $n);
	
	## GENERATES N RANDOM NUMBERS BETWEEN 0-NUMBER OF ENTRIES
	my %hash;
	for(my $i=0; $i<$n; $i++){
		my $flag = 0;	
		while($flag == 0){
			my $random_number = int(rand($range));
			if(exists $hash{$random_number}){
				$flag = 0;
			}else{
				$hash{$random_number} = $random_number;
				$flag = 1;
			}
		}
	}
	
	my $counter = 0;
	foreach(sort{$a <=>$b} keys %hash){
		$counter++;
	}
	#print STDERR "Numbers generated: ".$counter."\n";
		
	# Then print file(s).
	$counter = 0;
	if($infile){
		my $ref_fastq_db = Iterator::FastqDb->new($infile, {trimN=>0, trim3=>0}) or die("Unable to open Fastq file, $infile\n");
		while( my $curr = $ref_fastq_db->next_seq() ) {
			if(exists $hash{$counter}){
				print OUT $curr->output;
			}
			$counter++;
		}
	}elsif($infile_1 and $infile_2){
		my $ref_fastq_db_1 = Iterator::FastqDb->new($infile_1, {trimN=>0, trim3=>0}) or die("Unable to open Fastq file, $infile_1\n");
		my $ref_fastq_db_2 = Iterator::FastqDb->new($infile_2, {trimN=>0, trim3=>0}) or die("Unable to open Fastq file, $infile_2\n");
		while( my $curr1 = $ref_fastq_db_1->next_seq() ) {
			my $curr2 = $ref_fastq_db_2->next_seq();
			if(exists $hash{$counter}){
				print OUT_1 $curr1->output;
				print OUT_2 $curr2->output;
			}
			$counter++;
		}	
	}
}
close(OUT) if($outfile);
close(OUT_1) if($outfile_1);
close(OUT_2) if($outfile_2);

# Count number of lines in a file (it is assumed it is a fastq file) and divide by 4. Return the value.
# @input: one fastq file
# @output: a number
sub count_fastq_reads{
    my $file = $_[0];
    my $count = 0;

    if(-s $file and -e $file){
        $count = `cut -c -2 $file | wc -l`;
        chomp($count);
        $count =~ m/(\d+)/;
        $count = $1 / 4;
        #$count = commify($count);
    }
    #count = commify($count);

    return $count;
}

# Return commat formatted numbers
# $_[0] - a number
# @returns comma-formatted number
sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text
}

