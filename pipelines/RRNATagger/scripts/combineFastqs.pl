#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use lib "/global/homes/t/tremblay/build/Itagger/lib/perl5/JGI/Bio/Itagger/";
use Iterator::FastqDb;
use Iterator::FastaDb;
use File::Basename;
use Env qw/TMPDIR/;
use File::Temp;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
combine_multiple_fastq_libs.pl

PURPOSE:

INPUT:
Specify one:

--infile_fastq <string>       :    input fastq files, must be 2 or more. 
--infile_gzip <string>        :    input fastq files, must be 2 or more. 
--index <string>              :    barcodes file

OUTPUT:
--outfile_fastq <string>      :    new fastq outfile
--outfile_index <string>      :    barcodes in fasta format with:
   >new_name
      ACGTACGT
      .
      .
      .

--outfile_reference <string>  :    tabular index outfile with old(i.e.indexes) and "new" indexes.
   
   For instance, output will look like the following:
   <file name>\t<original_index_seq>\t<new_index_seq>\t<old_index_name>\t<new_index_name>\n
       run1.fastq	ACGTACGT	ACGTACGTACGT	12_sept_2012_site_A		12_sept_2012_site_A_run1.fastq
       .
       .
       .
       run2.fastq	ACGTACGA	ACGTACGAACGT	30_oct_2012_site_A		30_oct_2012_site_A_run2.fastq
NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, @infile_fastq, @infile_gzip, $outfile_index, $outfile_fastq, @index, $outfile_reference);
my $verbose = 0;

GetOptions(
    'infile_fastq=s' 		=> \@infile_fastq,
    'infile_gzip=s' 		=> \@infile_gzip,
	'index=s' 				=> \@index,
    'outfile_index=s' 		=> \$outfile_index,
    'outfile_fastq=s' 		=> \$outfile_fastq,
	'outfile_reference=s' 	=> \$outfile_reference,
    'verbose' 				=> \$verbose,
    'help' 					=> \$help
);
if ($help) { print $usage; exit; }

# TEMP FILES
my $tmpdir = File::Temp->newdir(
    "tmpDirCombineFastqsXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

# VALIDATION
die "Provide at least 2 --index\n" if(@index < 2);
die "Provide --outfile_fastq" unless($outfile_fastq);
die "Provide --outfile_index" unless($outfile_index);
die "Provide --outfile_reference" unless($outfile_reference);
if($infile_fastq[0]){
	die "Provide at least 2 --infile_fastq\n" if(@infile_fastq < 2);
	die "Specify equal number of --index and --infile_fastq args\n" if(@infile_fastq != @index);
}
if($infile_gzip[0]){
	die "Provide at least 2 --infile_gzip\n" if(@infile_gzip < 2);
	die "Specify equal number of --index and --infile_gzip args\n" if(@infile_gzip != @index);
}

foreach(@index){
	die("Barcode file does not exists... Wrong filename? : ".$_."\n") if(!-e $_ or !-s $_);
}


## MAIN
die "No tempdir specified, please specify a \$TMPDIR variable in your environment.\n" if(!$TMPDIR);

open(OUT_FASTQ, ">".$outfile_fastq) or die "Can't open output file ".$outfile_fastq."\n";
open(OUT_INDEX, ">".$outfile_index) or die "Can't open output file ".$outfile_index."\n";
open(OUT_REF, ">".$outfile_reference) or die "Can't open output file ".$outfile_reference."\n";


## GUNZIP .gz FILE
my $uncompressed_file_name;
my $infile;
my @fastqs;

## Here chech if files are gzips or fastqs
if($infile_gzip[0]){
	foreach(@infile_gzip){
		#print $tmpdir."\n";
		$uncompressed_file_name = basename($_);
		$uncompressed_file_name =~ s/\.fastq//;
		$uncompressed_file_name =~ s/\.gz//;
		$uncompressed_file_name =~ s/\.gzip//;
	    $uncompressed_file_name = $tmpdir."/".$uncompressed_file_name.".fastq";
    	print "Decompressing gzip to fastq";
		print $uncompressed_file_name."\n" if($verbose);
		print "Decompressing file:\t gunzip -c ".$_." > ".$uncompressed_file_name;
	    system("gunzip -c ".$_." > ".$uncompressed_file_name);
	    push(@fastqs, $uncompressed_file_name);
	}
}elsif($infile_fastq[0]){
	foreach(@infile_fastq){
		push(@fastqs, $_);
		print $_."\n" if($verbose);
	}
}

print "Finished unzipping files\n" if($verbose);

## GENERATE DUMMY INDEXES
my @new_barcodes = getBarcodes();

foreach(@new_barcodes){
	#print $_."\n" if($verbose);
}

my @test_index = @index;
my @test_fastqs = @fastqs;

## FIRST SAMPLE LENGTH OF BARCODES IN INDEX FILE AND FASTQ FILE
## We assume that withing a barcodes file, all barcodes are of the same length.
my $barcode_min_length = 1000;
foreach(@test_index){
	my $barcode_length;
	my $temp_barcode = new Iterator::FastaDb($_) or die("Unable to open Fastq file, $_\n");
	while( my $curr = $temp_barcode->next_seq() ){
		$barcode_length = length($curr->seq);
		last;	
	}
	$barcode_min_length = $barcode_length if($barcode_length < $barcode_min_length);
}
die "Barcode length can't be of 0. Please double-check your barcode fasta files...\n" if($barcode_min_length == 0);
print STDERR "BARCODE MIN LENGTH:\t ".$barcode_min_length."\n" if($verbose);

my $barcode_fastq_min_length = 1000;
foreach(@test_fastqs){
	my $barcode_length_fastq;
	my $temp_fastq = new Iterator::FastqDb($_) or die("Unable to open Fastq file, $_\n");
	while( my $curr = $temp_fastq->next_seq() ){
		print "CURR BARCODE:\t".$curr->barcode."\n";
		my $barcode = uc($curr->barcode());
		$barcode_length_fastq = length($barcode);
		last;
	}
	$barcode_fastq_min_length = $barcode_length_fastq if($barcode_length_fastq < $barcode_fastq_min_length);
	
	print STDERR "BARCODE LENGTH IN FASTQ:\t".$barcode_length_fastq."\n" if($verbose);
}

my $length;
if($barcode_min_length > $barcode_fastq_min_length){
	$length = $barcode_fastq_min_length;
}elsif($barcode_min_length < $barcode_fastq_min_length){
	$length = $barcode_min_length;
}else{
	$length = $barcode_min_length;
}
print STDERR "FINAL BARCODE LENGTH:\t".$length."\n";
die "Final barcode length must be > 0...\n" if($length < 1);

sleep(2);

# GENERATE 1 UNIQUE FILE WITH NEW UNIQUE AND LONGER BARCODES
my $reference_string = "#File_name\toriginal_barcode_seq\tnew_barcode_seq\told_barcode_name\tnew_barcode_name\n";

foreach my $fastq (@fastqs){
	my $lib_name = basename($fastq);
	my $index = shift(@index);
	
	print STDERR "Processing $fastq\t" if($verbose);
	print STDERR "Processing $index\n" if($verbose);
	
	#STORE ORIGINAL INDEXES IN A HASH;
	my @old_barcodes = ();
	my @old_barcodes_names = ();
	
	my $fasta_in = new Iterator::FastaDb($index);
	while( my $curr = $fasta_in-> next_seq() ){
		push(@old_barcodes, uc(substr($curr->seq, 0, $length)));
		push(@old_barcodes_names, substr($curr->header(), 1));
		print $curr->seq()."\n".substr($curr->header(),1)."\n" if($verbose);
	}
	
	my %hash = ();
	#ADD NEW INDEX VALUE TO ORIGINAL INDEX_VALUE
	foreach my $old_barcode (@old_barcodes){
		my $new_barcode = shift(@new_barcodes); #here new barcodes are permanently removed from @new_barcodes array.
		my $old_barcode_name = shift(@old_barcodes_names);
		print $old_barcode."\t".$new_barcode."\n" if($verbose);
		$hash{$old_barcode} = uc($new_barcode);
		$reference_string .= $lib_name."\t".$old_barcode."\t".$new_barcode."\t".$old_barcode_name."\t".$old_barcode_name."_".$lib_name."\n";
		print OUT_INDEX ">".$old_barcode_name."_".$lib_name."\n".$new_barcode."\n"; #Will write new index to file even if there are no matching sequence in file.
	}

	# DEFINE VARIANTS OF BARCODE (allow single base sequencing error)
	#my $barcodes = barcode_variants(@old_barcodes);	


	#LOOP THROUGH FASTQ FILE, REPLACE OLD BARCODE WITH NEW BARCODE FOR EACH SEQ AND WRITE NEW FILE.
	my $in = new Iterator::FastqDb($fastq) or die("Unable to open Fastq file, $fastq\n");
	while( my $curr = $in->next_seq() ){
		my $barcode = uc(substr($curr->barcode, 0, $length));
		my $pair = $curr->pair();
		my $base = $curr->base();
		
		#If barcode exists, replace it

		if(exists $hash{$barcode}){
		#if($mismatch_counter <= 1){
			print OUT_FASTQ "@".$base."#".$hash{$barcode}."/".$pair."\n".$curr->seq()."\n+\n".$curr->qual()."\n";
		}	
	}
}

print OUT_REF $reference_string;

close(OUT_REF);
close(OUT_INDEX);
close(OUT_FASTQ);

sub getBarcodes {

	my $i=0;
	my $j=0;
	my @barcodes;
	for(glob '{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}{A,C,G,T}'){
	    if($i % 30 == 0){
			$barcodes[$j] = $_;
			$j++;
		}else{
		
		}
		$i++;
	}
	
	## CHECK IF EACH BARCODES ARE UNIQUE.
	my %hash = ();
	foreach(@barcodes){
		$hash{$_} = 1;
	}
	
	#foreach(@barcodes){
	my %hash_check = ();
	while ( my ($key, $value) = each(%hash) ) {
		if(exists $hash_check{$key}){
			$hash_check{$key}++;
			die("Barcodes is in duplicate: ".$key."\n");
		}else{
			$hash_check{$key} = 1;
		}
	}
	
	my @barcodes_nodup;
	while ( my ($key, $value) = each(%hash_check) ) {
		if($value <= 1){
			#print "$key => $value\n" 
			push(@barcodes_nodup, $key);
		}
	}
	
	return @barcodes_nodup;	
	#exit;
}


sub reverse_complement_IUPAC {
        my $dna = shift;

    # reverse the DNA sequence
        my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}

## REMOVE TEMP FILES
sub END{
	system("rm ".$tmpdir." -rf");
}

