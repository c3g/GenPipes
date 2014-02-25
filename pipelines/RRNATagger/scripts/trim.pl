#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use Bio::Seq;
use Bio::Tools::IUPAC;
use File::Which;
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
Trim.pl

PURPOSE:

INPUT:
--infile <string>    :  Sequence file
--outfile <string>   :  Sequence file
--num_threads <int>  :  Number of threads

--meanq <int>        :  default=15
--winsize <int>      :  default=3
--minlen <int>       :  default=75
--maxn <int>         :  default=20 

OUTPUT:

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $num_threads, $meanq, $winsize, $minlen, $maxn);

my $verbose = 0;

GetOptions(
    'infile=s' 			=> \$infile,
	'outfile=s' 		=> \$outfile,
	'meanq=i' 			=> \$meanq,
	'winsize=i' 		=> \$winsize, 
	'minlen=i' 			=> \$minlen,
	'maxn=i' 			=> \$maxn,
	'num_threads=i' 	=> \$num_threads,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--outfile arg required\n") unless($outfile);
die("--num_threads arg required\n") unless($num_threads);

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerTrimXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

#my jgi_qc = which('itagger_illumina_read_qc.pl');
#if(!defined($jgi_qc)) {
#	die "itagger_illumina_read_qc.pl is not on the path...\n";
#}else{
#    print STDOUT "Using : ".$jgi_qc."\n" if($verbose);
#}

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
#my @out_fhs_failed;

for(my $i=0; $i < $num_threads; $i++){
	push(@out_fhs, $tmpdir."/"."_trimmed_".$i);
	#push(@out_fhs_failed, $outfile_failed."_".$i);
}

## Calculate each threads start posn
my $quarter = int( $size / $num_threads );
my @starts = map $quarter * $_, 0 .. ($num_threads - 1);
push @starts, $size;

## Start N threads and wait for them to end.
$_->join for map
{
    async( \&worker, $infile, @starts[ $_, $_ +1 ], $out_fhs[$_])
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

    $tmpdir = $tmpdir."/$$".$tid;
    mkdir $tmpdir or die("Can't make dir ".$tmpdir."\n");
    print $tmpdir."\n" if($verbose);

    open my $FASTQ, '<', $file or die $!;
    open my $FASTQ_TEMP, '>', $tmpdir."/fastq.temp" or die $!;

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
		
		if(defined ($pair) and defined ($barcode)){
			print $FASTQ_TEMP $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");
		}elsif(defined($pair)and !defined($barcode)){
			print $FASTQ_TEMP $base."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");
		}elsif(!defined($pair) and defined($barcode)){
			print $FASTQ_TEMP $base."#".$barcode."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");	
		}else{
			print $FASTQ_TEMP $base."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");
		}
		
    }
	close($FASTQ);
	close($FASTQ_TEMP);
	
	## Perform task here 	
	function($tmpdir."/fastq.temp", $outfile);	
}

sub function{
	my($infile, $unpaired_outfile) = @_;
	my $append = 0;

	# QC PARAMETERS
	my %qc_params = (
		'winsize' => $winsize,
		'meanq'   => $meanq,
		'minlen'  => $minlen,
		'maxn'    => $maxn,
		'low'     => 0.9,
	    'trim3'   => 0,
	    'trimN'   => 0
	);
	
	#my $unpaired = $unpaired_outfile ? new IO::File("$append$unpaired_outfile") : *STDOUT;
	open(OUT, ">".$unpaired_outfile);
	my %summary=();
	#while (@ARGV){
		#my $infile = shift @ARGV;
		my $db = new Iterator::FastqDb($infile, \%qc_params);
		my $read1;
		my $read2;
		while ($read2 = $db->next_seq){
			if ($read2->filtered){
				$read2 = undef;
			}elsif(!$read2->pair){
				if(length($read2->seq) >= $minlen){
		            print OUT $read2->output or die("NFS failure\n");
				}
				$read2 = undef;
			}elsif(defined($read1) and $read1->base eq $read2->base){
				# complete pair
				if(length($read1->seq) >= $minlen){	
					print OUT $read1->output, $read2->output or die("NFS failure\n");
				}
				$read1 = $read2 = undef;
			}
	
			if(defined($read1)){
				#$read1->unpair;
				if(length($read1->seq) >= $minlen){	
					print OUT $read1->output or die("NFS failure\n");
				}
				$read1 = undef;
			}
			$read1 = $read2;
		}
	    my $summary0 = $db->summary;
	    foreach my $key (keys %$summary0){
	        my $val = $summary0->{$key};
	        if( exists($summary{$key}) ){
	            $summary{$key} += $val;
	        }else{
				$summary{$key} = $val;
	        }
	    }
		close(OUT);
	#}
	#close($unpaired) if $unpaired_outfile;
}
exit;

## REMOVE TEMP FILES
sub END{
	local $?;
	system("rm ".$tmpdir." -rf");
}

