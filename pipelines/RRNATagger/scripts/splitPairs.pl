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

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
SplitPairs.pl

PURPOSE:
Parallel wrapper for itagger_separate_pairs.pl

INPUT:
--infile <string>    :  Sequence file
--outfile_1 <string> :  Sequence file
--outfile_2 <string> :  Sequence file
--num_threads <int>  :  Number of threads

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $outfile_1, $outfile_2, $num_threads);

my $verbose = 0;

GetOptions(
    'infile=s' 		=> \$infile,
	'outfile_1=s' 	=> \$outfile_1,
	'outfile_2=s' 	=> \$outfile_2,
	'num_threads=i' => \$num_threads,

    'verbose' 		=> \$verbose,
    'help' 			=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--outfile_1 arg required\n") unless($outfile_1);
die("--outfile_2 arg required\n") unless($outfile_2);
die("--num_threads arg required\n") unless($num_threads);

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerSplitPairsXXXXXXX",
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

for(my $i=0; $i < $num_threads; $i++){
	push(@out_fhs_1, $tmpdir."/_R1_".$i);
	push(@out_fhs_2, $tmpdir."/_R2_".$i);
}

## Calculate each threads start posn
my $quarter = int( $size / $num_threads );
my @starts = map $quarter * $_, 0 .. ($num_threads - 1);
push @starts, $size;

## Start N threads and wait for them to end.
$_->join for map
{
    async( \&worker, $infile, @starts[ $_, $_ +1 ], $out_fhs_1[$_], $out_fhs_2[$_] )
} 0 .. ($num_threads - 1);

## Concatenate file parts.
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
$cat .= "> ".$outfile_2;
print $cat."\n" if($verbose);
system($cat);
system("rm ".$_) foreach(@out_fhs_2);

print time - $t."\t".$num_threads."\n" if($verbose);

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
    my( $file, $start, $end, $outfile_1, $outfile_2) = @_;

	print $outfile_1."\n" if($verbose);
	print $outfile_2."\n" if($verbose);
    open my $FASTQ, '<', $file or die $!;
	open my $FASTQ_OUT_1, '>', $outfile_1 or die $!;
	open my $FASTQ_OUT_2, '>', $outfile_2 or die $!;
	#open my $FASTQ_OUT_FAILED, '>', $outfile_failed or die $!;

    ## If a no-zero start posns, find the start of the next full record.
	## Here the $FASTQ file handler will be modified in the findNextRecStart (implicit return).
    findNextRecStart( $FASTQ, $start, $tid ) if $start;
	
    ## process records until the end of this threads section.
    while( tell( $FASTQ ) < $end ) {
        my @lines = map scalar( <$FASTQ> ), 1 .. 4;
        chomp @lines;
	
		## Validate header
		my $validate = new Iterator::ValidateFastq($lines[0], $lines[3], $phred);
		my $base = "@".$validate->base;
		my $barcode = $validate->barcode;
		my $pair = $validate->pair;
		my $header = $lines[0];
		my $seq = $lines[1];
		#my $qual = $lines[3];
		my $qual = $validate->qual;	
	
		## Perform task here 
		die("No paring information is available in the fastq header\n") unless($pair);	
		function($header, $base, $barcode, $pair, $seq, $qual, $FASTQ_OUT_1, $FASTQ_OUT_2);	
    }
	close $FASTQ_OUT_1;
	close $FASTQ_OUT_2;
}

sub function{
	my($header, $base, $barcode, $pair, $seq, $qual, $OUT_1, $OUT_2) = @_;
	if ( !defined($pair) ) {
        die($header." does not having pairing information\n");
    } elsif ( $pair == 1 ) {
		if(defined $barcode){
        	print $OUT_1 $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n";
		}else{
        	print $OUT_1 $base."/".$pair."\n".$seq."\n+\n".$qual."\n";	
		}
    } elsif ($pair == 2) {
		if(defined $barcode){
        	print $OUT_2 $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n";
		}else{
        	print $OUT_2 $base."/".$pair."\n".$seq."\n+\n".$qual."\n";
		}
    } else {
        die($header." has invalid pair ID, ".$seq."\n");
    }
}

## REMOVE TEMP FILES
sub END{
	system("rm ".$tmpdir." -rf");
}
