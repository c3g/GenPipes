#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastqDb;

my $usage=<<'ENDHERE';
NAME:
Preprocess454.pl

PURPOSE:

INPUT:
--fasta <string>         :  input fasta files, can be 1 and more. 
--qual <string>          :  qual files, can be 1 and more
--index <string>         :  a tab format barcode sequence file. Ex: name\tsequence\n, can be 1 and more
--length <int>           :  minimum length of reads that will be in output.
--outdir <string>        :  output directory.

OUTPUT:
--outfile_fastq <string> : fastq outfile
--outfile_index <string> : index outfile with "new" indexes

NOTES:
Must have multiple --fasta and --qual and --index and there must be equal number of both and be in the 
same order. For instance jgi_Preprocess454.pl \
--fasta 1.fasta --fasta 2.fasta \
--qual 1.qual --qual 2.qual \
--index 1.primers --index 2.primers \

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, @fasta, $outfile_index, $outfile_fastq, @index, @qual, $outdir, $length);
my $verbose = 0;

GetOptions(
    'fasta=s' 			=> \@fasta,
	'qual=s' 			=> \@qual,
	'index=s' 			=> \@index,
	'length=i' 			=> \$length,
	'outdir=s' 			=> \$outdir,
    'outfile_index=s' 	=> \$outfile_index,
    'outfile_fastq=s' 	=> \$outfile_fastq,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATION
die "Provide at least 1 --qual\n" unless(@qual);
die "Provide at least 1 --fasta\n" unless(@fasta);
die "Provide at least 1 --index\n" unless(@index);
die "Provide --outdir\n" unless($outdir);
die "Provide --outfile_fastq" unless($outfile_fastq);
die "Provide --outfile_index" unless($outfile_index);
die "Specify equal number of --qual, --index and --fasta args\n" if(@qual != @fasta and @qual != @index and @index != @fasta);

## MAIN
open(OUT_FASTQ, ">".$outfile_fastq) or die "Cant open file ".$outfile_fastq."\n";
open(OUT_INDEX, ">".$outfile_index) or die "Cant open file ".$outfile_index."\n";
open(OUT_INDEX_TAB, ">".$outfile_index.".tab") or die "Cant open file ".$outfile_index.".tab\n";
open(OUT_INDEX_COR, ">".$outfile_index.".cor") or die "Cant open file ".$outfile_index.".cor\n";

my @fastq_with_barcodes = ();
my @original_indexes = @index;
# ====================================== #
# Convert Fastq+Qual to Fastqs========== #
# ====================================== #
my $j=1;
foreach my $fasta (@fasta){
	my $qual = shift(@qual);
	my $index = shift(@index);
	system("fasta_qual_to_fastq ".$fasta." ".$qual." ".$outdir."/".$j.".fastq -s");
	print $fasta."\n" if($verbose);
	print $qual."\n" if($verbose);
	print $index."\n" if($verbose);

	my %primers = ();
	
	open(INDEX, $index);
	while(<INDEX>){
	    chomp($_);
	    my @row = split(/\t/, $_);
	    $primers{substr($row[1],0,5)} = 0;
	}
	close(INDEX);
	
	my $infile = $outdir."/".$j.".fastq";
	open(OUT_WITH_BARCODES, ">".$outdir."/".$j."_with_barcodes.fastq");

	my $in_1 = new Iterator::FastqDb($infile) or die("Unable to open Fastq file, $infile\n");
	while( my $curr = $in_1->next_seq() ){
	
	    my $header = "";
	    my $seq = "";
    	my $qual = "";

	    my $begin = substr($curr->seq, 0, 8);
	
	    #Barcode filter
		my $matched = 0;
	    for my $key (keys %primers){
	        if($begin =~ m/$key/g){
	            my $pos = pos($begin);
				$pos = $pos + 15;
	            $header = $curr->header()."#".$key;
	            $seq = substr($curr->seq(),$pos);
	            $qual = substr($curr->qual(),$pos);
				$matched = 1;
				last;
	        }
	    }
	
		if($matched == 1){
    		#Length filter
		    if(defined($length)){
				if(length($seq) >= $length){
		        	print OUT_WITH_BARCODES $header."\n".substr($seq,0,$length)."\n+\n".substr($qual,0,$length)."\n";
	    		}#Else, dont print the sequence...
			}else{
		    	print OUT_WITH_BARCODES $header."\n".$seq."\n+\n".$qual."\n";
		    }
		}
	}
	push(@fastq_with_barcodes, $outdir."/".$j."_with_barcodes.fastq");
	close(OUT_WITH_BARCODES);
	$j++;
}
#At this point, good fastq file in proper format are generated.

# ====================================== #
# Generate dummy primers================ #
# ====================================== #
my @array = ();
my $i=0;
my $k=0;
for(glob '{A,C,G,T,T}{A,C,G,T,T}{A,C,G,T,T}{A,C,G,T,T}{A,C,G,T,T}{A,C,G,T,T}{A,C,G,T,T}'){
    if($i % 10 == 0){
        $array[$k] = $_;
        $k++;
    }else{

    }
    $i++;
}

## CHECK IF EACH BARCODES ARE UNIQUE.
my %hash = ();
foreach(@array){
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

my @dummy_primers;
while ( my ($key, $value) = each(%hash_check) ) {
   # print "$key => $value\n" if($value <= 1);
	push(@dummy_primers, $key);
}

# ====================================== #
# Generate 1 unique file with new unique # 
# and longer barcodes ================== #
# ====================================== #
my %id = ();
$i=0;
foreach my $file (@original_indexes){
	my %primers = ();
	my $fastq = shift(@fastq_with_barcodes);
	print $file."\n".$fastq."\n" if($verbose);
	
	open(INDEX, $file);
	while(<INDEX>){
		chomp($_);
		my @row = split(/\t/, $_);
		$primers{substr($row[1],0,5)} = $dummy_primers[$i];
		print OUT_INDEX ">".$row[0]."\n".$dummy_primers[$i]."\n";	

		$id{$dummy_primers[$i]} = substr($row[1],0,5)."\t".$row[0];
		$i++;
	}
	close(INDEX);
	
	#Write new file with new barcodes here
	my $in = new Iterator::FastqDb($fastq) or die("Unable to open Fastq file, $fastq\n");
	while( my $curr = $in->next_seq() ){
		my @header = split(/#/,$curr->header());
		
		if(exists $primers{$header[1]}){
			print OUT_FASTQ $header[0]."#".$primers{$header[1]}."/1\n".$curr->seq()."\n+\n".$curr->qual()."\n";
		}	
	}
	$fastq = "";
	$in = "";
}

## Print final references files.
for my $key (keys %id){
	print OUT_INDEX_COR $key."\t".$id{$key}."\n";
	my @row = split(/\t/, $id{$key});
	print OUT_INDEX_TAB $row[1]."\t".$key."\n";	
}
exit;

