#!/usr/bin/env perl

use strict;
use warnings;

use Symbol 'delete_package';
use Env qw/TMPDIR/;
use String::Approx 'aslice';
use List::Util qw(sum);
use threads;
use threads::shared;
use Getopt::Long;
use File::Temp;
use Parallel::ForkManager;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use Cache::FastMmap;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
ParallelRDP.pl

PURPOSE:
Parallel wrapper for RDP classifier.

INPUT:
--infile <string>                : fastq (or fasta) formatted file
--rdp_training_set <string>      : File path pointing to RDP training set.
--minWords <int>                 : Set the minimal number of words to use. New with RDP classifier
                                   2.5. Default = 120.
--fasta                          : If file is in fasta format. 
--num_threads <int>              : Number of threads.

OUTPUT:
--outfile <string>               : outfile containing rdp values.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $rdp_training_set, $minWords, $num_threads, $fasta);

my $verbose = 0;

GetOptions(
    'infile=s' 				=> \$infile,
	'rdp_training_set=s' 	=> \$rdp_training_set,
	'minWords=i'			=> \$minWords,
	'outfile=s' 			=> \$outfile,
	'num_threads=i' 		=> \$num_threads,
	'fasta'					=> \$fasta,
    'verbose' 				=> \$verbose,
    'help' 					=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--outfile arg required\n") unless($outfile);
die("--num_threads arg required\n") unless($num_threads);
die("--rdp_training_set arg required\n") unless($rdp_training_set);
$minWords = 120 unless($minWords);
die("--minWords must be <= 200 \n") if($minWords > 200);

## TEMP DIR
my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerRDPXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

die "RDP rRNAClassifier.properties file does not exist or is empty...\nFILE: ".$rdp_training_set."\n"  if((!-e $rdp_training_set) or (!-s $rdp_training_set));

## CHECK RDP CLASSIFIER PATH
my $rdp_classifier = `which rdp_classifier.jar`; die if(!defined $rdp_classifier); chomp $rdp_classifier;


####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################

my $t = time;
my $phred;
if($fasta){
	my $count = 0;
	my $ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fastq file, $infile\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		$count++;
	}

	my $seq_per_query = $count/$num_threads;
	$seq_per_query = int($seq_per_query) + 1;
	print $seq_per_query."\n" if($verbose);
	
	my %hash = ();
	my @array = ();
	my $counter = 0;
	my $file_counter = 1;
	my $curr_file;
	
	$curr_file = $tmpdir."/seq".$file_counter.".fasta";

	open(CURR, ">".$curr_file);
	#$hash{$out_file} = $file_counter if($out_file);
	push(@array, $curr_file);


	$ref_fasta_db = Iterator::FastaDb->new($infile) or die("Unable to open Fasta file, $infile\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
	
		print CURR $curr->header()."\n".$curr->seq()."\n";
	
		if($counter == $seq_per_query){
			close(CURR);
			$counter = 0;
			$file_counter++;
			$curr_file = $tmpdir."/seq".$file_counter.".fasta";
			#print $curr_file."\n" if($verbose);
			push(@array, $curr_file);
			open(CURR, ">".$curr_file);
	}
		$counter++;
	}
	#push(@array, $curr_file);
	close(CURR);
	
	# Start Fork manager
	my $pm = new Parallel::ForkManager($num_threads);
	my $index = 0;
	my $Cache = Cache::FastMmap->new();
	$Cache->set(0, "");	

	foreach(@array){
		$index++;
		my $pid = $pm->start($index) and next;

		## LAUNCH THE RDP CLASSIFIER
		my $cmd = "java -Xmx1g -jar ".$rdp_classifier;
		$cmd .= " -q ".$_;
		$cmd .= " -o ".$tmpdir."/OUT_".$index;
		$cmd .= " -t ".$rdp_training_set;
		$cmd .= " --minWords ".$minWords;
		`$cmd  2>&1`;
		$? != 0 ? die "command failed: $cmd\n" : print STDOUT "RDP classifier (inside ParallelRDP.pl script) successfuly executed\n" if($verbose);
		print $cmd."\n" if($verbose);	

		$Cache->set($index-1, $tmpdir."/OUT_".$index);	

		# Terminate ForkManager
		$pm->finish($index);
	}
	$pm->wait_all_children;

	# concat all RDP sheets.
	my $cmd = "cat";
	for(my $i=0;$i<$num_threads;$i++){
		$cmd .= " ".$Cache->get($i);
	}
	$cmd .= " > ".$outfile;
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "Concataneting RDP sheets done\n" if($verbose);	

	
# If fastq.
}else{
	## Unload Cache::FastMmap module because not compatible with threads.
	delete_package 'Cache::FastMmap';

	## Check if qual is phred 33 or phred 64.
	$phred = Iterator::Utils->new($infile);
	
	## Grab the size of the file
	my $size =  -s $infile or die "$! : $infile";
	#my $num_threads = $num_threads or die "$! : $num_threads";
	
	my @out_fhs;
	
	for(my $i=0; $i < $num_threads; $i++){
		push(@out_fhs, $outfile."_".$i);
	}
	
	## Calculate each threads start posn
	my $quarter = int( $size / $num_threads );
	my @starts = map $quarter * $_, 0 .. ($num_threads - 1);
	push @starts, $size;
	
	## Start N threads and wait for them to end.
	$_->join for map
	{
	    async( \&worker, $infile, @starts[ $_, $_ +1 ], $out_fhs[$_] )
	} 0 .. ($num_threads - 1);
	
	## Concatenate file parts.
	my $cat = "cat ";
	foreach(@out_fhs){
		$cat .= $_. " ";
	}
	$cat .= "> ".$outfile;
	print $cat."\n" if($verbose);
	system($cat);
	system("rm ".$_) foreach(@out_fhs);
	
}
cleanup();

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
    my( $file, $start, $end, $outfile) = @_;

	print $outfile."\n" if($verbose);
    open my $FASTQ, '<', $file or die $!;
	open my $FASTA_TEMP, '>', $tmpdir."/temp_fasta_".$tid.".fasta" or die $!;
	#open my $FASTQ_TEMP, '>', $tmpdir."/temp_fastq_".$tid.".fastq" or die $!;
	#print $tmpdir."/temp_fastq_".$tid.".fastq\n";

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
		#my $pair = $validate->pair;
		my $header = $lines[0];
		my $seq = $lines[1];
		#my $qual = $lines[3];
		my $qual = $validate->qual;
		
		## Perform task here 
		print $FASTA_TEMP ">".$base."#".$barcode."\n".$seq."\n";
		#print $FASTQ_TEMP $lines[0]."\n".$lines[1]."\n+\n".$lines[3]."\n";
    }
	close $FASTA_TEMP;
	#close $FASTQ_TEMP;
	
	## LAUNCH THE RDP CLASSIFIER
	my $cmd = "java -Xmx1g -jar ".$rdp_classifier;
	$cmd .= " -q ".$tmpdir."/temp_fasta_".$tid.".fasta";
	$cmd .= " -o ".$outfile;
	$cmd .= " -t ".$rdp_training_set;
	$cmd .= " --minWords ".$minWords;
	`$cmd  2>&1`;
	print STDERR $cmd."\n" if($verbose);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "RDP classifier (inside ParallelRDP.pl script) successfuly executed\n" if($verbose);	
}

sub function{
	my($header, $base, $barcode, $pair, $seq, $qual, $OUT_1, $OUT_2) = @_;
	# DO SOMETHING COOL HERE...
}

## REMOVE TEMP FILES
sub cleanup{
    system("rm ".$tmpdir." -rf");
}
