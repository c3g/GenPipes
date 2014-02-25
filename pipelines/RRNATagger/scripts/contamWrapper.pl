#!/usr/bin/env perl

use strict;
use warnings;

#use lib "/global/homes/t/tremblay/build/JGI/JGI/perl/Itagger/lib";
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
contamWrapper.pl

PURPOSE:
Parallel wrapper for the duk program.

INPUT:
--infile <string>            : Sequence file
--db <string>                : Reference db in fasta format
--num_threads <int>          : Number of threads

--kmer <int>                 : Default = 21
--step <int>                 : Default = 1
--cutoff <int>               : Default = 1
	
OUTPUT:
--outfile_matched <string>   : Matched sequences
--outfile_unmatched <string> : Unmatched sequences
--log <string>               : Log file

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $debug, $infile, $matched, $unmatched, $db, $num_threads, $log, $kmer, $step, $cutoff); 

my $verbose = 0;

GetOptions(
    'infile=s' 				=> \$infile,
	'num_threads=i' 		=> \$num_threads,
	'log=s'					=> \$log,
	'db=s'					=> \$db,
    'outfile_matched=s'		=> \$matched,
    'outfile_unmatched=s'	=> \$unmatched,
	'kmer=i'				=> \$kmer,
	'step=i'				=> \$step,
	'cutoff=i'				=> \$cutoff,
    'verbose' 				=> \$verbose,
    'help' 					=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n")            unless($infile);
die("--db arg required\n")                unless($db);
die("--log arg required\n")               unless($log);
die("--outfile_matched arg required\n")   unless($matched);
die("--outfile_unmatched arg required\n") unless($unmatched);
die("--num_threads arg required\n")       unless($num_threads);
die("--infile file is empty or does not exists! (Typed wrong filename?)\n") if((!-e $infile) and (!-s $infile));
####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################

$kmer = 21 unless($kmer);
$step = 1 unless($step);
$cutoff = 1 unless($cutoff);

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerContamXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

my $t = time;
share($db);

## Check if qual is phred 33 or phred 64.
my $phred = Iterator::Utils->new($infile);

## Grab the size of the file
my $size =  -s $infile or die "$! : $infile";
#my $num_threads = $num_threads or die "$! : $num_threads";

my @out_fhs;
my @out_fhs_failed;
my @out_log;

for(my $i=0; $i < $num_threads; $i++){
	push(@out_fhs, $matched."_".$i);
	push(@out_fhs_failed, $unmatched."_".$i);
	push(@out_log, $log."_".$i);
}

## Calculate each threads start posn
my $quarter = int( $size / $num_threads );
my @starts = map $quarter * $_, 0 .. ($num_threads - 1);
push @starts, $size;

## Start 4 threads and wait for them to end.
$_->join for map
{
    async( \&worker, $infile, @starts[ $_, $_ +1 ], $out_fhs[$_], $out_fhs_failed[$_], $out_log[$_], $tmpdir)
} 0 .. ($num_threads - 1);

## Concatenate file parts.
my $cat = "cat ";
foreach(@out_fhs){
	$cat .= $_. " ";
}
$cat .= "> ".$matched;
print $cat."\n" if($verbose);
system($cat);
system("rm ".$_) foreach(@out_fhs);

## Concat failed
$cat = "cat ";
foreach(@out_fhs_failed){
	$cat .= $_. " ";
}
$cat .= "> ".$unmatched;
print $cat."\n" if($verbose);
system($cat);
system("rm ".$_) foreach(@out_fhs_failed);

## Concat logs
$cat = "cat ";
foreach(@out_log){
	$cat .= $_. " ";
}
$cat .= "> ".$log;
print $cat."\n" if($verbose);
system($cat);
system("rm ".$_) foreach(@out_log);

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
    #twarn "[$tid] $skipLines";

    ## Seek bck to that calculated start posn.
    seek $fh, $start, 0;
	#print "Start/End after found:\t".$startFound."\n";

    ## Then skip forward th calculate dnumber of lines.
    scalar <$fh> for 1 .. $skipLines;
	#twarn "[$tid] ", tell $fh;
		
    return;
}

sub worker {
    my $tid = threads->tid;
    ## the name of the file, the byte offsets for the thread
    ## to start and end processing at
    my( $file, $start, $end, $outfile, $outfile_failed, $log, $tmpdir) = @_;

	print $outfile."\n" if($verbose);
    open my $FASTQ, '<', $file or die $!;
	open my $FASTQ_TEMP, '>', $tmpdir."/duk_temp_".$tid or die $!;

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

		if ( defined($pair) && defined($barcode) ) {
        	print $FASTQ_TEMP $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n";
	    } elsif ( defined($pair) && !defined($barcode) ) {
        	print $FASTQ_TEMP $base."/".$pair."\n".$seq."\n+\n".$qual."\n";	
	    }else{
        	print $FASTQ_TEMP $base."#".$barcode."\n".$seq."\n+\n".$qual."\n";
		} 			
    }
	close $FASTQ_TEMP;

	## Execute duk
	my $cmd = "duk";
	$cmd .= " -o ".$log;
	$cmd .= " -n ".$outfile_failed;
	$cmd .= " -m ".$outfile;
	$cmd .= " -k ".$kmer;
	$cmd .= " -s ".$step;
	$cmd .= " -c ".$cutoff;
	$cmd .= " ".$db;
	$cmd .= " ".$tmpdir."/duk_temp_".$tid;
	system($cmd);
	die "command failed: $!\n" if($? != 0);
}

## REMOVE TEMP FILES
sub END{
	local $?;
	system("rm ".$tmpdir." -rf");
}
