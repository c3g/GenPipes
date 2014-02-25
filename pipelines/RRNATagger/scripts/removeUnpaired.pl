#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use String::Approx 'aslice';
use List::Util qw(sum);
use File::Temp;
use threads;
use threads::shared;
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use Devel::Size qw(size total_size);

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
RemoveUnpaired.pl

PURPOSE:
Remove unpaired reads from a paired end reads fastq file.
Paired end reads of a fastq file might get disrupted
after filtering for contaminants. This scripts assumes
that reads are in a collated paired end reads order
For instance:
@A/1
@A/2
@B/2
@C/1
@C/2
@D/1
Will only keep valid pairs having both /1 and /2.

INPUT:
--infile <string>         :  Sequence file
--outfile_1 <string>      :  Sequence file
--outfile_2 <string>      :  Sequence file
--outfile_paired <string> :  Sequence file
--num_threads <int>       :  Number of threads

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $outfile_1, $outfile_2, $outfile_paired, $num_threads);

my $verbose = 0;

GetOptions(
    'infile=s' 			=> \$infile,
	'outfile_1=s' 		=> \$outfile_1,
	'outfile_2=s' 		=> \$outfile_2,
	'outfile_paired=s' 	=> \$outfile_paired,
	'num_threads=i' 	=> \$num_threads,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") 			unless($infile);
die("--outfile_1 arg required\n") 		unless($outfile_1);
die("--outfile_2 arg required\n") 		unless($outfile_2);
die("--outfile_paired arg required\n") 	unless($outfile_paired);
die("--num_threads arg required\n") 	unless($num_threads);

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDirRmUnpairedXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################
my $t = time;

## Check if qual is phred 33 or phred 64.
my $phred = Iterator::Utils->new($infile);

## Grab the size of the file
my $size =  -s $infile or die "$! : $infile";
#my $num_threads = $num_threads or die "$! : $num_threads";

my @out_fhs_1;
my @out_fhs_2;
my @out_fhs_paired;
my %master_hash_1 :shared;
my %master_hash_2 :shared;

## Prepare output files and array of hashes
for(my $i=0; $i < $num_threads; $i++){
	push(@out_fhs_1, $tmpdir."/_pair1_".$i);
	push(@out_fhs_2, $tmpdir."/_pair2_".$i);
	push(@out_fhs_paired, $tmpdir."/_paired_".$i);
}

## Calculate each threads start posn
my $quarter = int( $size / $num_threads );
my @starts = map $quarter * $_, 0 .. ($num_threads - 1);
push @starts, $size;

# Start N threads and wait for them to end.
$_->join for map
{
    async( 
		\&worker, 
		$infile, 
		@starts[ $_, $_ +1 ], 
		$out_fhs_1[$_], 
		$out_fhs_2[$_], 
		$out_fhs_paired[$_], 
		\%master_hash_1,
		\%master_hash_2 
	);
} 0 .. ($num_threads - 1);

# Concatenate file parts.
my $cat = "cat ";
foreach(@out_fhs_1){
	$cat .= $_. " ";
}
$cat .= "> ".$outfile_1;
print $cat."\n" if($verbose);
system($cat);
system("rm ".$_) foreach(@out_fhs_1);

$cat = "cat ";
foreach(@out_fhs_2){
	$cat .= $_. " ";
}
$cat .= "> ".$outfile_2.".temp";
print $cat."\n" if($verbose);
system($cat);
system("rm ".$_) foreach(@out_fhs_2);

$cat = "cat ";
foreach(@out_fhs_paired){
	$cat .= $_. " ";
}
$cat .= "> ".$outfile_paired;
print $cat."\n" if($verbose);
system($cat);
system("rm ".$_) foreach(@out_fhs_paired);

# look into hashes and see if there are mate pairs.
open my $FASTQ_OUT_1, '>>', $outfile_1 or die $!;
open my $FASTQ_OUT_2, '>', $outfile_2 or die $!;
open my $FASTQ_OUT_PAIRED, '>>', $outfile_paired or die $!;

my %hash_flag_2;

while ( my ($key, $value) = each(%master_hash_1) ) {
	#print "$key => $value\n";
	if(exists $master_hash_2{$key}){
		print $FASTQ_OUT_PAIRED $value;
		print $FASTQ_OUT_PAIRED $master_hash_2{$key};
		delete $master_hash_1{$key};
		$hash_flag_2{$key} = $master_hash_2{$key};
		delete $master_hash_2{$key};	
	}
}

## Populate unpaired /1 fastq file
while ( my ($key, $value) = each(%master_hash_1) ) {
	print $FASTQ_OUT_1 $value;
}

## Last step is to remove /2 reads that might be in duplicate. Skip if no pairs /2
if( (-e $outfile_2.".temp") and (-s $outfile_2.".temp")){
	my $db = Iterator::FastqDb->new($outfile_2.".temp", {trimN=>0, trim3=>0}) or die("Unable to open Fastq file, ".$outfile_2.".temp\n");
	while( my $curr = $db->next_seq() ) {
		my $base = "@".$curr->base;
		if(exists $hash_flag_2{$base}){
			## Don't print seqs flagged in the hash		
		}else{
			print $FASTQ_OUT_2 $curr->output;
		}
	}
}

close $FASTQ_OUT_1;
close $FASTQ_OUT_2;
close $FASTQ_OUT_PAIRED;

system("rm ".$outfile_2.".temp");

print "Time: ".(time - $t)."\tThreads:".$num_threads."\n" if($verbose);

exit;

#################
## SUBROUTINES ##
#################

## threadsafe output routines.
$|++; ## Does work without this!
my $semStdout :shared;
sub tprint{ lock $semStdout; print @_; }
my $semStderr :shared;
sub twarn{ lock $semStderr; print STDERR @_; }

sub findNextRecStart {
    ## filehandle, calculated start byte, thread id (for tracing
    my( $fh, $start, $tid ) = @_;
	#twarn "[$tid] $start";

    ## seek to the start byte -1; Just incase the calculated posn hits bang on
    seek $fh, $start-1, 0;

    ## Read a buffer full; we'd be really unluck to not find a record start in 4k
    ## But you could increase this to say 64k if it fails.
    read( $fh, my $buffer, 4096 );

    ## Search for a full record that doesn't have @/+ as the first char in the 2nd line
    $buffer =~ m[\n(\@)(?:[^@+][^\n]+\n){2}\+] or die "Couldn't locate record start";

    ## Remember the offset into the buffer where we found it.
    my $startFound = $-[1];

    ## Now count the lines between the start of the buffer and that po int.
	my @ASCII = unpack("C*", substr( $buffer, 0, $startFound )); #Added fix when what is separating header's '@' char from previous line is equal to \n only (ASCII 10).
    my $previousLines = substr( $buffer, 0, $startFound ) =~ tr[\n][\n];
    if( @ASCII == 1 && $ASCII[0] == 10 ){ #Added fix: so when only a line feed is separating header's '@' char from previous line, put $previous line to 0 instead of 1.
        $previousLines = 0;
    }elsif(@ASCII > 1 && $ASCII[0] == 10){
		$previousLines = $previousLines - 1;
	}

    ## And calulate our way back to the first full record after the calculated start posn.
    my $skipLines = ( $previousLines - 1) % 4 +1;
    #twarn "[$tid] $skipLines";

    ## Seek bck to that calculated start posn.
    seek $fh, $start, 0;

    ## Then skip forward th calculate dnumber of lines.
    scalar <$fh> for 1 .. $skipLines;
	#twarn "[$tid] ", tell $fh;
    return;
}

sub worker {
    my $tid = threads->tid;
    ## the name of the file, the byte offsets for the thread
    ## to start and end processing at
    my( $file, $start, $end, $outfile_1, $outfile_2, $outfile_paired, $hash_1, $hash_2) = @_;

    open my $FASTQ, '<', $file or die $!;
	open my $FASTQ_OUT_1, '>', $outfile_1 or die $!;
	open my $FASTQ_OUT_2, '>', $outfile_2 or die $!;
	open my $FASTQ_OUT_PAIRED, '>', $outfile_paired or die $!;

    ## If a no-zero start posns, find the start of the next full record.
	## Here the $FASTQ file handler will be modified in the findNextRecStart (implicit return).
    findNextRecStart( $FASTQ, $start, $tid ) if $start;
	
	my $pair_1=0;
	my $pair_2=0;
	my $counter=0;
	my $loop_counter=0;
	my $seq_1;
	my $seq_2;
	my $seq_1_base;
	my $seq_2_base;

	my $base;
	my $barcode;
	my $pair;
	my $header;
	my $seq;
	my $qual;
	my $output;

    ## process records until the end of this threads section.
	while( tell( $FASTQ ) < $end ) {
		my @lines = map scalar( <$FASTQ> ), 1 .. 4;
		chomp @lines;

		## Validate header
		my $validate = new Iterator::ValidateFastq($lines[0], $lines[3], $phred);
		$base = "@".$validate->base;	
		$barcode = $validate->barcode;	
		$pair = $validate->pair;
		$header = $lines[0];
		$seq = $lines[1];
		#$qual = $lines[3];
		$qual = $validate->qual;
	
		## Perform task here
		unless($pair){
			print STDERR "================\n";
			print STDERR $header."\n".$seq."\n+\n".$qual."\n";
			print STDERR "================\n";
			die("No paring information is available in the fastq header\n");
		}
		$output = $header."\n".$seq."\n+\n".$qual."\n";
		#function($header, $base, $barcode, $pair, $seq, $qual, $FASTQ_OUT_1, $FASTQ_OUT_2, $FASTQ_OUT_PAIRED);

		## If chunk of file starts with a /2 seq, perhaps the /1 pair is in the
		## Previous file chunk.
		$master_hash_2{$base} = $output if(($loop_counter == 0) && ($pair == 2) );
		$loop_counter++;

	    if ( !defined($pair) ) {
	        die($seq->header." does not having pairing information\n");
			next;
	    } elsif($pair == 1 ) {
			$pair_1++;
			$seq_1 = $output;
			$seq_1_base = $base;
			$counter++;
	    } elsif($pair == 2) {
			$pair_2++;
			$seq_2 = $output;
			$seq_2_base = $base;
			$counter++;
	    } else {
	        die($header." has invalid pair ID, ".$pair."\n");
	    }
	
		#validate	
		if($pair_1 == 1 && $pair_2 == 1 && $counter == 2 && $seq_1_base eq $seq_2_base){ #GOOD PAIRS!
			print $FASTQ_OUT_PAIRED $seq_1.$seq_2;
			$pair_1=0;
			$pair_2=0;
			$counter=0;
			$seq_1="";
			$seq_2="";
			$seq_1_base=0;
			$seq_2_base=0;
		}elsif($pair_1 == 1 && $pair_2 == 1 && $counter == 2 && $seq_1_base ne $seq_2_base){ #GOOD PAIRS, but base id does not correspond
			print $FASTQ_OUT_1 $seq_1;
			print $FASTQ_OUT_2 $seq_2;
			$pair_1=0;
			$pair_2=0;
			$counter=0;
			$seq_1="";
			$seq_2="";
			$seq_1_base=0;
			$seq_2_base=0;
		}elsif($pair_1 == 1 && $pair_2 == 0 && $counter == 2){ #Not sure this can happen...
			print $FASTQ_OUT_1 $seq_1;
			$pair_1=0;
			$pair_2=0;
			$counter=0;
			$seq_1="";
			$seq_1_base=0;
			$seq_2_base=0;
			next;	
		}elsif($pair_1 == 0 && $pair_2 == 1 && $counter == 2){ #Not sure this can happen...
			print $FASTQ_OUT_2 $seq_2;
			$pair_1=0;
			$pair_2=0;
			$counter=0;
			$seq_2="";
			$seq_1_base=0;
			$seq_2_base=0;
			next;	
		}elsif($pair_1 == 0 && $pair_2 == 1 && $counter == 1){ #Problem, because pair 2 comes before pair 1. Means that there is no pair one.
			print $FASTQ_OUT_2 $seq_2;
			$pair_1=0;
			$pair_2=0;
			$counter=0;
			$seq_2="";
			#$seq_1_base=0;
			$seq_2_base=0;
			next;
		#}elsif($pair_1 == 1 && $pair_2 == 0 && $counter == 1){
		#	print $FASTQ_OUT_1 $seq_1;
		#	$pair_1=0;
		#	$pair_2=0;
		#	$counter=0;
		#	$seq_1="";
		#	$seq_1_obj=0;
		#	#$seq_2_obj=0;
		#	next;		
		}elsif($pair_1 == 2 && $counter == 2){
			print $FASTQ_OUT_1 $seq_1;
			$pair_1=1;
			$pair_2=0;
			$counter=1;
			$seq_2="";
			#$seq_1_base=0;
			$seq_2_base=0;
			next;	
		}else{
			#DEBUGGING
			if($verbose){
				print "PAIR 1: ".$pair_1."\n";
				print "PAIR 2: ".$pair_2."\n";
				print "COUNTER: ".$counter."\n";
				print $seq_1."\n";
				print $seq_2."\n";
			}
		}
	}
	
	# Put the last sequence in corresponding hash.
	if($pair == 1){
		$master_hash_1{$base} = $output;
	}
	close $FASTQ;
	close $FASTQ_OUT_1;
	close $FASTQ_OUT_2;
	close $FASTQ_OUT_PAIRED;
}

## REMOVE TEMP FILES
sub END{
	local $?;
	system("rm ".$tmpdir." -rf");
}
