#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use Bio::Seq;
use Bio::Tools::IUPAC;
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
use Data::Dumper;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
itagsQC.pl

PURPOSE:
Parallel wrapper for Itags QC

INPUT:
--infile <string>          : Sequence file
--num_threads <int>        : Number of threads
--primer_5_prime <string>  : Fasta file containing ONE primer sequence
--primer_3_prime <string>  : Fasta file containing ONE primer sequence
--length_5_prime <int>     : Length of the 5' region to search for primer
--length_3_prime <int>     : Length of the 3' region to search for primer 
--primer_mismatch <int>    : % of mismatch allowed when searching for primer.
                             Default: 15.
--cut_first <int>          : Cut the first <int> bases
--cut_last <int>           : Cut last <int> bases
--qscore_1 <int>           : Average Qscore.
--N <int>                  : Max number of N tolerated
--lq_threshold <int>       : Low quality threshold
--qscore_2 <int>           : Number of tolerated bases below low quality 
                             threshold
--max_length <int>         : Sequences having more than <int> bases length
                             will be rejected.
--min_length <int>         : Sequences having less than <int> bases length
                             will be rejected. 
--qual <int>               : Phred scale. Either 33 or 64.
--noqc                     : Set if no QC is to be done. For instance if primers
                             are to be removed, but no QC to be done.
--remove_1st_N             : Trim first base of read in 5' region if it is a N
                             (Optional).
--num_threads <int>        : Number of threads to use.

OUTPUT:
--outfile_failed <string>  : File containing reads that failed to pass
                             QC.
--reject_unmatched         : Set this flag if sequences not matching
                             primer seq. To use in combination with
                             either --primer_5_prime or --primer_3_prime

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $debug, $infile, $outfile, $outfile_failed, $num_threads,
    $primer_5_prime, $primer_3_prime, $length_5_prime, $length_3_prime,
    $primer_mismatch, $cut_first, $cut_last, $qscore_1, $N, $lq_threshold,
    $qscore_2, $max_length, $min_length, $qual, $noqc, $reject_unmatched,
	$remove_N
);
my $verbose = 0;

GetOptions(
    'infile=s' 			=> \$infile,
	'outfile=s' 		=> \$outfile,
	'num_threads=i' 	=> \$num_threads,
    'infile=s' 			=> \$infile,
    'primer_5_prime=s' 	=> \$primer_5_prime,
    'primer_3_prime=s' 	=> \$primer_3_prime,
    'length_5_prime=i' 	=> \$length_5_prime,
    'length_3_prime=i' 	=> \$length_3_prime,
    'primer_mismatch=i' => \$primer_mismatch,
	'remove_1st_N'		=> \$remove_N,
    'cut_last=i' 		=> \$cut_last,
    'cut_first=i' 		=> \$cut_first,
    'qscore_1=i' 		=> \$qscore_1,
    'N=i' 				=> \$N,
    'lq_threshold=i' 	=> \$lq_threshold,
    'qscore_2=i' 		=> \$qscore_2,
    'max_length=i' 		=> \$max_length,
    'min_length=i' 		=> \$min_length,
    'qual=i' 			=> \$qual,
    'outfile_failed=s' 	=> \$outfile_failed,
    'reject_unmatched' 	=> \$reject_unmatched,
    'noqc' 				=> \$noqc,
    'debug' 			=> \$debug,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--infile $infile might be empty or wrong file path?\n") if((!-e $infile) and (!-s $infile));
die("--outfile arg required\n") unless($outfile);
die("--num_threads arg required\n") unless($num_threads);
if($primer_mismatch){
    die("--primer_mismatch int value between 0 and 50 required\n") unless($primer_mismatch > 0 and $primer_mismatch < 51);
}
if($max_length){
    die("--max_length Please enter a value <= 1000\n") if($max_length > 1000);
}
die("--outfile outfile required\n") unless $outfile;

if( ($primer_5_prime and !$length_5_prime) or (!$primer_5_prime and $length_5_prime) ){
    die "--primer_5_prime AND --length_5_prime must be provided\n"
}else{
    $primer_mismatch = 15 unless $primer_mismatch;
}
if( ($primer_3_prime and !$length_3_prime) or (!$primer_3_prime and $length_3_prime) ){
    die "--primer_3_prime AND --length_3_prime must be provided\n"
}else{
    $primer_mismatch = 15 unless $primer_mismatch;
}
$primer_mismatch = 15 unless $primer_mismatch;

my $left_primer = 0;
my $right_primer = 0;
if($primer_mismatch && $primer_5_prime && $length_5_prime){
    $left_primer = 1;
}
if($primer_mismatch && $primer_3_prime && $length_3_prime){
    $right_primer = 1;
}
if($noqc){
    $noqc = 1;
}else{
    $noqc = 0;
}
if($noqc == 0){
    $qscore_1 = 30 unless $qscore_1;
    $N = 3 unless $N;
    $lq_threshold = 3 unless $lq_threshold;
    $qscore_2 = 15 unless $qscore_2;
    die "please enter either 33 or 64 for --qual arg\n" unless $qual == 33 or $qual == 64;
}

$cut_first = 0 unless($cut_first);
$cut_last = 0 unless($cut_last);

## Temp files
my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerItagsQcXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

## PREPARE PRIMERS
my @left_primer;
my @right_primer;

## GENERATE PRIMER SEQS
if($primer_5_prime){
    my $left_primer_seq;
    my $ref_fasta_db = new Iterator::FastaDb($primer_5_prime) or die("Unable to open Fasta file, $primer_5_prime\n");
    my $i=0;
    while( my $ref_seq = $ref_fasta_db->next_seq() ) {
        $left_primer_seq = uc($ref_seq->seq());
        $i++;
    }
    die "Please provide only one sequence in fasta file\n" if($i > 1);
    @left_primer = generate_primer_combination($left_primer_seq);
}

if($primer_3_prime){
    my $right_primer_seq;
    my $ref_fasta_db = new Iterator::FastaDb($primer_3_prime) or die("Unable to open Fasta file, $primer_3_prime\n");
    my $i=0;
    while( my $ref_seq = $ref_fasta_db->next_seq() ) {
        $right_primer_seq = uc($ref_seq->seq());
        $right_primer_seq =~ tr/ACGT/TGCA/; #Just convert to complement since the reads will be reversed
        $i++;
    }
    die "Please provide only one sequence in fasta file\n" if($i > 1);
    @right_primer = generate_primer_combination($right_primer_seq);
}
#my %qscore_hash :shared;
#setQscoreHash(\%qscore_hash);
#print Dumper(\%qscore_hash);

####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################
my $t = time;

## Check if qual is phred 33 or phred 64.
my $phred = Iterator::Utils->new($infile);

## Grab the size of the file
my $size =  -s $infile or die "$! : $infile";
#my $num_threads = $num_threads or die "$! : $num_threads";

my @out_fhs;
my @out_fhs_failed;

for(my $i=0; $i < $num_threads; $i++){
	push(@out_fhs, $tmpdir."/_passed_".$i);
	push(@out_fhs_failed, $tmpdir."/_failed_".$i);
}

## Calculate each threads start posn
my $quarter = int( $size / $num_threads );
my @starts = map $quarter * $_, 0 .. ($num_threads - 1);
push @starts, $size;

## Start 4 threads and wait for them to end.
$_->join for map
{
    async( \&worker, $infile, @starts[ $_, $_ +1 ], $out_fhs[$_], $out_fhs_failed[$_] )
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

$cat = "cat ";
foreach(@out_fhs_failed){
	$cat .= $_. " ";
}
$cat .= "> ".$outfile_failed;
print $cat."\n" if($verbose);
system($cat);
system("rm ".$_) foreach(@out_fhs_failed);

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

    ## Seek bck to that calculated start posn.
    seek $fh, $start, 0;

    ## Then skip forward th calculate dnumber of lines.
    scalar <$fh> for 1 .. $skipLines;
		
    return;
}

sub worker {
    ## the name of the file, the byte offsets for the thread
    ## to start and end processing at
    my( $file, $start, $end, $outfile, $outfile_failed) = @_;
    my $tid = threads->tid;
	print "Current thread:\t".$tid."\n" if($verbose);

	print $outfile."\n" if($verbose);
    open my $FASTQ, '<', $file or die $!;
	open my $FASTQ_OUT, '>', $outfile or die $!;
	open my $FASTQ_OUT_FAILED, '>', $outfile_failed or die $!;

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
		
		filter($lines[0],$lines[1],$lines[3], $FASTQ_OUT, $FASTQ_OUT_FAILED);
		#print $lines[0]."\n".$lines[1]."\n+\n".$lines[3]."\n";
			
    }
	close $FASTQ_OUT;
	close $FASTQ_OUT_FAILED;
}

## FILTER SEQUENCE
## input : sequence_header_string, dna_string, quality_string, out_filehandler, out_failed_filehandler.
## output: null
sub filter{
	my($header, $newseq, $newqual, $OUT, $OUT_F) = @_;
	my %qscore_hash;
	setQscoreHash(\%qscore_hash);

    #REMOVE LEFT PRIMER
    if($left_primer == 1){
		my $string = substr($newseq, 0, $length_5_prime);
        my $matched = 0;
        my $lowest_length = 1000; #Arbitrary unrealistic high value.

        foreach(@left_primer){
            my @index = aslice($_, ["i ".$primer_mismatch."%"] , $string);
            if(defined $index[0][0]){
    	            my $length_from_0 = $index[0][0] + $index[0][1];
                $lowest_length = $length_from_0 if($length_from_0 < $lowest_length);
    	            $matched = 1;
        	    }
        }

		if($matched == 1 and $lowest_length < 1000){
            $newseq = substr($newseq, $lowest_length);
            $newqual = substr($newqual,$lowest_length);
        }else{
            if($reject_unmatched){
                print $OUT_F $header."\n".$newseq."\n+\n".$newqual."\n";
                next; #Discard read.
            }
        }
    }

    #REMOVE RIGHT PRIMER
    if($right_primer == 1){
        #Reverse both seq/qual
        $newseq = reverse($newseq);
        $newqual = reverse($newqual);
        my $string = substr($newseq, 0, $length_3_prime);
        my $matched = 0;
        my $lowest_length = 1000; #Arbitrary unrealistic high value.
 	
		foreach(@right_primer){
			
			print "REV Primer:\t".$_."\n" if($debug);

            my @index = aslice($_, ["i ".$primer_mismatch."%"] , $string);

            if(defined $index[0][0]){
				my $length_from_0 = $index[0][0] + $index[0][1];
                $lowest_length = $length_from_0 if($length_from_0 < $lowest_length);
                $matched = 1;
            }
        }

        if($matched == 1 and $lowest_length < 1000){
            $newseq = substr($newseq, $lowest_length);
   	        $newqual = substr($newqual,$lowest_length);
        }else{
            if($reject_unmatched){
                print $OUT_F $header."\n".$newseq."\n+\n".$newqual."\n";
                next; #discard
            }
        }
        print "Right primer Header:\t".$header."\n" if $debug;
		#print substr($newseq, 0, 20)."\n";

        $newseq = reverse($newseq);
        $newqual = reverse($newqual);
    }
    
	#Cut sequences before going into quality threshold.
	if($cut_first or $cut_last){
		my $length = (length($newseq)) - $cut_last - $cut_first;
		die "First nucleotides to cut must be higher or equal than/to read length.\n" if($cut_first >= length($newseq));
		$newseq  = substr($newseq, $cut_first, $length);
		$newqual = substr($newqual, $cut_first, $length);
	}
	
	#Remove 1st base if it is a N (typically seen in 5'reads).
	if($remove_N){
		if(substr($newseq, 0, 1) eq "N"){
			$newseq = substr($newseq, 1);	
			$newqual = substr($newqual, 1);	
		}
		
		# Also remove last base if last base is a N.	
		if(substr($newseq, -1) eq "N"){
			$newseq = substr($newseq, 0, -1);	
			$newqual = substr($newqual, 0, -1);	
		}
	
	}

	#QUALITY CONTROL (REMOVE READS WITH TOO MUCH Ns AND TOO MUCH NT HAVING LOW QUAL SCORE)
    if($noqc == 0){
        my @qual = unpack("C*", $newqual);
		my $prob_sum = 0;
		foreach(@qual){
			print $_ - $qual."\t" if($debug);
			$prob_sum = $prob_sum + $qscore_hash{$_ - $qual};
			die "Q score value does not exists in reference hash...".$_ - $qual."\n" if(!exists $qscore_hash{$_ - $qual});
		}
		
		die "Sequence has length of 0 bases. It you supplied the cut_last or cut_last paramters, it is possible that it made the sequence too short. Try to cut shorter...\n" if(@qual <= 0);			

		my $average = $prob_sum/@qual;
		print "\nAverage:\t".$average."\n" if($debug);
		
        my $Q_count = 0;
        foreach(@qual){
            $Q_count++ if ($_ - $qual) < $qscore_2;
        }
        print STDERR "Q_count:\t".$Q_count."\n" if $debug;

        my @seq = split(//, $newseq);
        my $N_count = 0;
        foreach(@seq){
            $N_count++ if uc($_) eq "N";
        }
        print $N_count."\n" if $debug;

        if( 	 $average > $qscore_hash{$qscore_1},
             and $N_count < $N,
             and $Q_count < $lq_threshold ){

            print "Average Q score is:\t".$average."\n" if $debug;

            if($max_length){
                #print "MAXLENGTH LOOP\n";
                if( length($newseq) <= $max_length ){
                     print STDERR $header."max_length passed\n" if $debug;
                }else{
                    next;
                }
			}elsif($min_length){
                if( length($newseq) >= $min_length ){
					print $OUT $header."\n".$newseq."\n+\n".$newqual."\n";
					print STDERR $header."min_length passed\n" if $debug;
                }else{
                    #next;
                }	
            }else{
				print $OUT $header."\n".$newseq."\n+\n".$newqual."\n";
                print STDERR $header."no min_length and no max_length\n" if $debug;
            }

        #Reject read if it does not encounter filtering parameters
        }else{
			print $OUT_F $header."\n".$newseq."\n+\n".$newqual."\n";
        }

	
	# If no QC is to be done.
    }else{
        if($cut_first or $cut_last){
			my $length = (length($newseq)) - $cut_last - $cut_first;
			print $OUT $header."\n".substr($newseq, $cut_first, $length)."\n+\n".substr($newqual, $cut_first, $length)."\n";
        }else{
            print $OUT $header."\n".$newseq."\n+\n".$newqual."\n";
        }
    }
}

## GENERATE ALL POSSIBLE PRIMER COMBINATION
## input: one dna sequence that can contain ambiguous bases.
## output: an array of dna sequences with no ambiguous bases.
sub generate_primer_combination{
    my $primer_seq = $_[0];
    my @array = ();

    #For some obscur reason, Bio::Seq would always accept sequence as provided...
    $primer_seq =~ m/(\S+)/;
    $primer_seq = $1;

    my $ambiseq = Bio::Seq->new(-seq => $primer_seq, -alphabet => 'dna');
    my $stream  = Bio::Tools::IUPAC->new(-seq => $ambiseq);

    while (my $uniqueseq = $stream->next_seq()) {
        push(@array, $uniqueseq->seq());
    }
    return @array;
}

## SET QSCORE HASH
## Will populate a hash ref passed in paramter.
## input: a hash reference
## output: null;
sub setQscoreHash{
    my($hash_ref) = @_;
#	$hash_ref->{0}=5.00;
#	$hash_ref->{1}=20.567;
#	$hash_ref->{2}=36.904;
#	$hash_ref->{3}=49.881;
#	$hash_ref->{4}=60.189;
#	$hash_ref->{5}=68.377;
#	$hash_ref->{6}=74.881;
#	$hash_ref->{7}=80.047;
#	$hash_ref->{8}=84.151;
#	$hash_ref->{9}=87.411;
#	$hash_ref->{10}=90.000;
#	$hash_ref->{11}=92.057;
#	$hash_ref->{12}=93.690;
#	$hash_ref->{13}=94.988;
#	$hash_ref->{14}=96.019;
#	$hash_ref->{15}=96.838;
#	$hash_ref->{16}=97.488;
#	$hash_ref->{17}=98.005;
#	$hash_ref->{18}=98.415;
#	$hash_ref->{19}=98.741;
#	$hash_ref->{20}=99.000;
#	$hash_ref->{21}=99.206;
#	$hash_ref->{22}=99.369;
#	$hash_ref->{23}=99.499;
#	$hash_ref->{24}=99.602;
#	$hash_ref->{25}=99.684;
#	$hash_ref->{26}=99.749;
#	$hash_ref->{27}=99.800;
#	$hash_ref->{28}=99.842;
#	$hash_ref->{29}=99.874;
#	$hash_ref->{30}=99.900;
#	$hash_ref->{31}=99.921;
#	$hash_ref->{32}=99.937;
#	$hash_ref->{33}=99.950;
#	$hash_ref->{34}=99.960;
#	$hash_ref->{35}=99.968;
#	$hash_ref->{36}=99.975;
#	$hash_ref->{37}=99.980;
#	$hash_ref->{38}=99.984;
#	$hash_ref->{39}=99.987;
#	$hash_ref->{40}=99.990;
#	$hash_ref->{41}=99.992;

	$hash_ref->{0}=1;
	$hash_ref->{1}=1;
	$hash_ref->{2}=2;
	$hash_ref->{3}=3;
	$hash_ref->{4}=4;
	$hash_ref->{5}=5;
	$hash_ref->{6}=6;
	$hash_ref->{7}=7;
	$hash_ref->{8}=8;
	$hash_ref->{9}=9;
	$hash_ref->{10}=10;
	$hash_ref->{11}=11;
	$hash_ref->{12}=12;
	$hash_ref->{13}=13;
	$hash_ref->{14}=14;
	$hash_ref->{15}=15;
	$hash_ref->{16}=16;
	$hash_ref->{17}=17;
	$hash_ref->{18}=18;
	$hash_ref->{19}=19;
	$hash_ref->{20}=20;
	$hash_ref->{21}=21;
	$hash_ref->{22}=22;
	$hash_ref->{23}=23;
	$hash_ref->{24}=24;
	$hash_ref->{25}=25;
	$hash_ref->{26}=26;
	$hash_ref->{27}=27;
	$hash_ref->{28}=28;
	$hash_ref->{29}=29;
	$hash_ref->{30}=30;
	$hash_ref->{31}=31;
	$hash_ref->{32}=32;
	$hash_ref->{33}=33;
	$hash_ref->{34}=34;
	$hash_ref->{35}=35;
	$hash_ref->{36}=36;
	$hash_ref->{37}=37;
	$hash_ref->{38}=38;
	$hash_ref->{39}=39;
	$hash_ref->{40}=40;
	$hash_ref->{41}=41;

}

## REMOVE TEMP FILES
sub END{
	local $?;
	system("rm ".$tmpdir." -rf");
}
