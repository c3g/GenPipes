#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use File::Temp;
use threads;
use threads::shared;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
Dereplicate.pl

PURPOSE:

INPUT:
--fasta <string>   : Sequence file (fasta format).
OR
--fastq <string>   : Sequence file (fastq format).
--sort             : Set if sequences are to be sorted and 
                     not dereplicated
--minsize <int>    : sequences below of minsize <int> will be
                     discarded.
				
OUTPUT:
--outfile <string> : Sequence file 100% clustered.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $fasta, $fastq, $outfile, $sort, $minsize, $num_threads);
my $verbose = 0;

GetOptions(
    'fasta=s' 		=> \$fasta,
	'fastq=s'		=> \$fastq,
	'sort'			=> \$sort,
	'minsize=i' 	=> \$minsize,
	'outfile=s'		=> \$outfile,
    'verbose' 		=> \$verbose,
	'num_threads=i' => \$num_threads,
    'help' 			=> \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--fasta OR --fastq required.\n" if(!defined($fastq) && !defined($fasta));
die "--outfile required.\n" unless($outfile);
if($fastq){
	die("--fastq does not exists! Typed wrong filename?\n") if((!-e $fastq) and (!-s $fastq));
}elsif($fasta){
	die("--fasta does not exists! Typed wrong filename?\n") if((!-e $fasta) and (!-s $fasta));
}

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerDereplicateXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0
);

$num_threads = 1 unless($num_threads);

####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################
my %hash = ();
my %hash_barcodes = ();
my %hash_headers = ();

my $t = time;

## Check if qual is phred 33 or phred 64.
my $phred = Iterator::Utils->new($fastq) if($fastq);

## If fastq file, split in multiple fasta before.
if($fastq){
	## Grab the size of the file
	my $size =  -s $fastq or die "$! : $fastq";
	#my $num_threads = $num_threads or die "$! : $num_threads";

	my @out_fhs;
	my @out_fhs_failed;
	
	for(my $i=0; $i < $num_threads; $i++){
		push(@out_fhs, $tmpdir."/file_".$i);
	}
	
	## Calculate each threads start posn
	my $quarter = int( $size / $num_threads );
	my @starts = map $quarter * $_, 0 .. ($num_threads - 1);
	push @starts, $size;
	
	## Start 4 threads and wait for them to end. Will write n_threads fasta files.
	$_->join for map
	{
	    async( \&worker, $fastq, @starts[ $_, $_ +1 ], $out_fhs[$_] )
	} 0 .. ($num_threads - 1);
	
	## Loop through fasta file parts and populate hash.
	constructHashFromFasta(@out_fhs);

}else{ ## Single threaded only.
	constructHashFromFasta($fasta);
}

sub constructHashFromFasta{

	my @files = @_;
	my $j=0;
	
	foreach my $file (@files){
		## FASTA
		my $ref_fasta_db = Iterator::FastaDb->new($file) or die("Unable to open Fasta file, $file\n");

		if($sort){# assumes that no 100% identical seqs are in dataset.
			while( my $curr = $ref_fasta_db->next_seq() ) {
				$curr->header =~ m/size=(\d+)/;
				my $size = $1;
				my $header = $curr->header;
				$header =~ s/size=\d+//;
	
				$hash{$curr->seq} = $size;
				$hash_headers{$curr->seq} = $header;
				
				print "processing seq ".$j."\n" if((($j % 100000) == 0) && $verbose );
				$j++;
			}
		}else{
			while( my $curr = $ref_fasta_db->next_seq() ) {
				if(exists $hash{$curr->seq}){
					$hash{$curr->seq}++;
					if(exists $hash_barcodes{$curr->seq}{$curr->barcode} ){
						$hash_barcodes{$curr->seq}{$curr->barcode}++;
					}else{
						$hash_barcodes{$curr->seq}{$curr->barcode}=1;
					}
				}else{
					$hash{$curr->seq} = 1;
					$hash_barcodes{$curr->seq}{$curr->barcode}=1;
				}
			
				print "processing seq ".$j."\n" if((($j % 100000) == 0) && $verbose );
				$j++;
			}
		}
	
	}
}

my $counter = 1;
open(OUT, ">".$outfile) or die "Can't open file ".$outfile."\n";

if($sort){ # Sort by size.
	$minsize = 0 unless($minsize);
	foreach my $key (sort {$hash{$b} <=> $hash{$a} } keys %hash){
		my $line = ">".$hash_headers{$key}.";size=".$hash{$key}."\n".$key."\n";
		$line =~ s/;;/;/;
		$line =~ s/>+/>/;
		print OUT $line if($hash{$key} >= $minsize);
	}
}else{
	foreach my $key (sort {$hash{$b} <=> $hash{$a} } keys %hash){
		print OUT ">".$counter.";";
		for my $barcode (keys %{ $hash_barcodes{$key} }){	
			print OUT "#".$barcode."=".$hash_barcodes{$key}{$barcode}.";";
		}
		print OUT "size=".$hash{$key}."\n";
		print OUT $key."\n";
		$counter++;
	}
}
close(OUT);

print time - $t."\t".$num_threads."\n" if($verbose);

exit;

#################
## SUBROUTINES ##
#################

## threadsafe output routines.
#$|++; ## Does work without this!
#my $semStdout :shared;
#sub tprint{ lock $semStdout; print @_; }
#my $semStderr :shared;
#sub twarn{ lock $semStderr; print STDERR @_; }

sub findNextRecStart {
    ## filehandle, calculated start byte, thread id (for tracing
    my( $fh, $start, $tid ) = @_;
	#twarn "[$tid] $start";

    ## seek to the start byte -1; Just incase the calculated posn hits bang on
    seek $fh, $start-1, 0;

    ## Read a buffer full; we'd be really unluck to not find a record start in 4k
    ## But you could increase this to say 64k if it fails.
    read( $fh, my $buffer, 8192 );

    ## Search for a full record that doesn't have @/+ as the first char in the 2nd line
    $buffer =~ m[\n(\@)(?:[^@+][^\n]+\n){2}\+] or die "Couldn't locate record start";

    ## Remember the offset into the buffer where we found it.
    my $startFound = $-[1];

	my @ASCII = unpack("C*", substr( $buffer, 0, $startFound )); #Added fix when what is separating header's '@' char from previous line is equal to \n only (ASCII 10).
    ## Now count the lines between the start of the buffer and that po int.
    my $previousLines = substr( $buffer, 0, $startFound ) =~ tr[\n][\n];
	if( @ASCII == 1 && $ASCII[0] == 10 ){ #Added fix: so when only a line feed is separating header's '@' char from previous line, put $previous line to 0 instead of 1.
		$previousLines = 0;
	}elsif(@ASCII > 1 && $ASCII[0] == 10){
		$previousLines = $previousLines - 1;
	}

    ## And calulate our way back to the first full record after the calculated start posn.
    my $skipLines = ( $previousLines - 1) % 4 +1;

    ## Seek bck to that calculated start posn.
    seek $fh, $start, 0;

    ## Then skip forward th calculate dnumber of lines.
    scalar <$fh> for 1 .. $skipLines;
		
    return;
}

sub worker {
    ## the name of the file, the byte offsets for the thread
    ## to start and end processing at
    my( $file, $start, $end, $outfile) = @_;
    my $tid = threads->tid;
	print "Current thread:\t".$tid."\n" if($verbose);

	print $outfile."\n" if($verbose);
    open my $FASTQ, '<', $file or die $!;
	open my $FASTA_OUT, '>', $outfile or die $!;

    ## If a no-zero start posns, find the start of the next full record.
	## Here the $FASTQ file handler will be modified in the findNextRecStart (implicit return).
    findNextRecStart( $FASTQ, $start, $tid ) if $start;

    ## process records until the end of this threads section.
    while( tell( $FASTQ ) < $end ) {
        my @lines = map scalar( <$FASTQ> ), 1 .. 4;
        chomp @lines;

		## Validate header
		my $validate = new Iterator::ValidateFastq($lines[0], $lines[3], $phred);
		#my $base = "@".$validate->base;	
		#my $barcode = $validate->barcode;	
		#my $pair = $validate->pair;
		my $header = $lines[0];
		my $seq = $lines[1];
		#my $qual = $lines[3];
		my $qual = $validate->qual;
		
		$header =~ s/@//;	
		print $FASTA_OUT ">".$header."\n".$seq."\n";
    }	
	close $FASTA_OUT;
}

print "Dereplicating seqs...\n" if($verbose);

print $minsize."\n" if($minsize and $verbose);

## REMOVE TEMP FILES
sub END{
	local $?;
	system("rm ".$tmpdir." -rf");
}

