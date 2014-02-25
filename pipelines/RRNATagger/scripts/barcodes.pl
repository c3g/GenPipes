#!/usr/bin/env perl

use strict;
no strict 'refs';
use warnings;

use Env qw/TMPDIR/;
use File::Which;
use String::Approx 'aslice';
use List::Util qw(sum);
use threads;
use threads::shared;
use Getopt::Long;
use Iterator::FastaDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use Data::Dumper;
use FileCache;
use File::Basename;
use File::Temp;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
Barcodes.pl

PURPOSE:

INPUT:
--infile <string>   : Sequence file.
--num_threads <int> : Number of threads.
--barcodes <string> : barcodes in fasta format.

OUTPUT:
--log <string>      : Log file having barcodes count.
--outfile <string>  : Sequence file in fastq having all good barcodes.
OR
--outdir <string>   : DIR where there will be one fastq file per barcode.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $outfile, $log, $outdir, $barcodes, $num_threads);

my $verbose = 0;

GetOptions(
    'infile=s' 			=> \$infile,
	'outfile=s' 		=> \$outfile,
	'outdir=s'			=> \$outdir,
	'barcodes=s'		=> \$barcodes,
	'log=s'				=> \$log,
	'num_threads=i' 	=> \$num_threads,
    'verbose' 			=> \$verbose,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
die("--infile arg required\n") unless($infile);
die("--infile might be empty or wrong file path?\n") if((!-e $infile) and (!-s $infile));
die("--outfile OR --outdir arg required\n") unless($outfile or $outdir);
die("--num_threads arg required\n") unless($num_threads);
die("Please specify a \$TMPDIR environment variable in your environment.\n") unless($TMPDIR);
die("Please specify a --log <string> arg/file.\n") unless($log);

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerBarcodesXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

my $outfile_failed;
if($outdir){
	$outfile_failed = $outdir."/UNKOWN.fastq";
}elsif($outfile){
	my($filename, $directories) = fileparse($outfile);
	$outfile_failed = $directories."/UNKNOWN_barcodes.fastq";
}else{
	die "Something went wrong...\n";
}

####################################
## MAIN - PREPARE SEQ REPARTITION ##
####################################
my $t = time;

## Check if qual is phred 33 or phred 64
my $phred = Iterator::Utils->new($infile);

## Grab the size of the file
my $size =  -s $infile or die "$! : $infile".". Cannot calculate size of file.\n";

my @out_fhs;
my @out_fhs_failed;
my @out_log;
my @total :shared;

my $outfile_base = basename($outfile) if($outfile);

for(my $i=0; $i < $num_threads; $i++){
	push(@out_fhs, $tmpdir."/".$outfile_base."_".$i) if($outfile);
	push(@out_log, $tmpdir."/barcodes_log_".$i);
	push(@total, 0);
	push(@out_fhs_failed, $tmpdir."/UNKNOWN_".$i);
}

## Calculate each threads start posn
my $quarter = int( $size / $num_threads );
my @starts = map $quarter * $_, 0 .. ($num_threads - 1);
push @starts, $size;

## Start N threads and wait for them to end.
$_->join for map
{
    async( \&worker, $infile, $barcodes, @starts[ $_, $_ +1 ], $out_fhs[$_], $out_fhs_failed[$_], $tmpdir, $out_log[$_], \@total)
} 0 .. ($num_threads - 1);

if($outfile){
	## Concatenate file parts.
	my $cat = "cat ";
	foreach(@out_fhs){
		$cat .= $_. " ";
	}
	$cat .= "> ".$outfile;
	print $cat."\n" if($verbose);
	system($cat);
	system("rm ".$_) foreach(@out_fhs);
	## Process each index file separately if $outdir.
}elsif($outdir){
	#print STDOUT "Concatenating\n";
	my $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file, $barcodes\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		my $header = $curr->header;
		$header =~ s/>//;
		$header =~ s/\//-/; # / will interfere in the moment of writing into file.
		$header =~ s/\s//g; # / will interfere in the moment of writing into file.
			
		my $cat_string = "cat ";
		my $found = "false";
		for(my $tid=1; $tid<=$num_threads; $tid++){
			next if(!-e $tmpdir."/".$header."_".$tid);
    		$cat_string .= $tmpdir."/".$header."_".$tid." ";
		}
		$cat_string .= " > ".$outdir."/".$header;
		
		#checkpoint to find if some barcodes files does not exists.
		if($cat_string =~ m/cat\s+> .*/){

		}else{
			#print "cat string:\t".$cat_string."\n";
			system($cat_string);
		
			# And then remove intermediate files.
			for(my $tid=1; $tid<=$num_threads; $tid++){
    			system("rm ".$tmpdir."/".$header."_".$tid);
			}
		}
	}
}

## Concatenate file parts for failed barcodes.
if($outdir){
	my $cat = "cat ";
	foreach(@out_fhs_failed){
		$cat .= $_. " ";
	}
	$cat .= "> ".$outfile_failed;
	print $cat."\n" if($verbose);
	system($cat);
	system("rm ".$_) foreach(@out_fhs_failed);
}

## Concatenate log files.
my %hash=();
my $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file, $barcodes\n");
while( my $curr = $ref_fasta_db->next_seq() ) {
	my $header = $curr->header;
	$header =~ s/>//;
	$header =~ s/\//-/; # / will interfere in the moment of writing into file.
	
	die "Barcode ".$curr->seq." is present in duplicate in barcode fasta file\n" if(exists $hash{$curr->seq});
	$hash{$curr->seq} = [$header, 0];
}

# Loop through log files
my $total = 0;
for(my $i=0; $i < $num_threads; $i++){
	#print "TOTAL:\t".$total[$i]."\n";
	$total = $total +  $total[$i];
	open(IN, "<".$tmpdir."/barcodes_log_".$i) or die "Can't open file ".$tmpdir."/barcodes_log_".$i;
	while(<IN>){
		chomp;
		my @row = split(/\t/, $_);
		##HEADER\tSEQUENCE\tCOUNT\n
		if(exists $hash{$row[1]}){
			$hash{$row[1]}[1] = $hash{$row[1]}[1] + $row[2];
		}else{
			##...
		}
	}
	close(IN);
}

# Print values in final log file.
open(LOG, ">".$log);
print LOG "#name\tsequence\tcount\tperc\n";
while ( my ($key, $value) = each(%hash) ) {
	my $perc = sprintf '%.2f', $hash{$key}[1]/$total * 100;
	print LOG $hash{$key}[0]."\t".$key."\t".$hash{$key}[1]."\t".$perc."%\n";
}
close(LOG);

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
	#print $tid."\n";
    ## the name of the file, the byte offsets for the thread
    ## to start and end processing at
    my( $file, $barcodes, $start, $end, $outfile, $outfile_failed, $outdir, $log, $total) = @_;

	## First, put barcodes into a hash.
	my %hash = ();
	my $ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file, $barcodes\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		my $header = $curr->header;
		$header =~ s/>//;
		$header =~ s/\//-/; # / will interfere in the moment of writing into file.
		
		die "Barcode ".$curr->seq." is present in duplicate in barcode fasta file\n" if(exists $hash{$curr->seq});
		$hash{$curr->seq} = [$header, 0];
		
		# If outdir, open file handlers for writing and push them into hash structure.
		if($outdir){
    		my $FASTQ_BARCODE = $outdir."/".$header."_".$tid;
			push @{$hash{$curr->seq}}, $FASTQ_BARCODE;
			
			# Open it once to initialize the file and overwrite potential old file.
    		open my $fh, '>', $FASTQ_BARCODE or die "Can't open ".$FASTQ_BARCODE." for writing.\n";
			close($fh);
		}
	}

	#print Dumper(\%hash);
    
	open my $FASTQ, '<', $file or die $!;
	my $FASTQ_OUT;
	if($outfile){
    	open $FASTQ_OUT, '>', $outfile or die $!;
	}
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
		my $qual = $validate->qual;
		#my $qual = $lines[3];
		
		## Perform barcode binning here.
		if(defined ($pair) and defined ($barcode)){
			if(exists $hash{$barcode}){
				$hash{$barcode}[1]++;
				if($outfile){
					print $FASTQ_OUT $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");		
				}elsif($outdir){
    				open my $fh, '>>', $hash{$barcode}->[2] or die $!;
					print $fh $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");		
					close($fh);
						
				}else{
					die "Somehting went wrong...\n";
				}
			}else{
				print $FASTQ_OUT_FAILED $base."#".$barcode."/".$pair."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");		
			}

		}elsif(defined($pair)and !defined($barcode)){
			die "Fastq sequence file does not contain barcodes\n";

		}elsif(!defined($pair) and defined($barcode)){
			if(exists $hash{$barcode}){
				$hash{$barcode}[1]++;
				if($outfile){
					print $FASTQ_OUT $base."#".$barcode."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");	
				}elsif($outdir){
    				open my $fh, '>>', $hash{$barcode}->[2] or die $!;
					print  $fh $base."#".$barcode."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");	
					close($fh);
	
				}else{
					die "Something went wrong...\n";
				}
			}else{
				print $FASTQ_OUT_FAILED $base."#".$barcode."\n".$seq."\n+\n".$qual."\n" or die("NFS failure\n");	
			}

		}else{
			die "Fastq sequence file does not contain barcodes and pair info\n";
		}	
    	
	}
	close($FASTQ);
	close($FASTQ_OUT) if($outfile);
	close($FASTQ_OUT_FAILED);
	
	open(LOG, ">".$log);
	while ( my ($key, $value) = each(%hash) ) {
		close($hash{$key}->[2]) if($outdir);
		##HEADER\tSEQUENCE\tCOUNT\n
		print LOG $hash{$key}[0]."\t".$key."\t".$hash{$key}[1]."\n";
		@$total[$tid-1] = @$total[$tid-1] +  $hash{$key}[1];
    }
	#print "TOTAL:\t".@$total[$tid-1]."\n";
	close(LOG);
	#print Dumper(\%hash);
}

## REMOVE TEMP FILES
sub END{
	local $?;
	system("rm ".$tmpdir." -rf");
}

