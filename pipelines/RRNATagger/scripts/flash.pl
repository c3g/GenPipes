#!/usr/bin/env perl

use strict;
use warnings;

use Env qw/TMPDIR/;
use List::Util qw(sum);
use Getopt::Long;
use File::Which;
use File::Temp;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Iterator::ValidateFastq;
use Iterator::Utils;
use Statistics::Descriptive;
#use SampleSheet;
use File::Basename;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
flash.pl

PURPOSE:
Wrapper for the FLASH software.

INPUT:
--infile_1 <string>            : Reads 1 (fastq)
--infile_2 <string>            : Reads 2 (fastq)
OR
--indir <string>               : If data to be assembled is from our GQIC internal
                                 libraries. The wrapper will look for the DIR 
                                 <INDIR>/raw_reads/
--sampleSheet <string>         : GQIC formatted sample sheet.

--num_threads <int>            : Number of threads
--n <int>                      : Number of reads you want to subsample to get
                                 paired-end assembly stats.
--computeQscores               : If quality stats plots are to be generated.
--prefix                       : Optional if suplying --infile_1 and --infile_2
                                 and you want to define your own prefix.

FLASH specific options
--m <int>                      : Min length default=15
--M <int>                      : Max length, default=280
--x <float>                    : Max mismatch ratio default=0.35
--p <int>                      : Quality scale (phred33 or phred64) default=33
--s <float>                    : standard deviation of assembled reads. default=5
--r <int>                      : read length default=150
--f <float>                    : length of assembled fragment default=150
	
OUTPUT:
--outdir <string>              : Output directory

NOTES:

BUGS/LIMITATIONS:
The flash binary has to be in your $PATH.

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile_1, $infile_2, $num_threads, $outdir, $start_at, $n, $sampleSheet, $indir, $computeQscores, $prefix,
	$m, $M, $x, $p, $s, $r, $f);

my $verbose = 0;

GetOptions(
    'infile_1=s' 		=> \$infile_1,
    'infile_2=s' 		=> \$infile_2,
	'indir=s'			=> \$indir,
	'sampleSheet=s'		=> \$sampleSheet,
	'outdir=s'			=> \$outdir,
	'prefix=s'			=> \$prefix,
	'num_threads=i' 	=> \$num_threads,
	'n=i'				=> \$n,
	'computeQscores'	=> \$computeQscores,
	'start_at=i'		=> \$start_at,

	'm=i' => \$m,
	'M=i' => \$M,
	'x=f' => \$x,
	'p=i' => \$p,
	's=f' => \$s,
	'r=i' => \$r,
	'f=f' => \$f,

    'verbose' => \$verbose,
    'help' => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE 
my $flash = which('flash');
if(!defined($flash)) {
    die("Exiting... The flash binary was not found on your \$PATH.\n");
}else{
    print STDERR "Using : ".$flash."\n" if($verbose);
}
die("--outdir arg missing...\n") unless($outdir);
if($indir){
	die("--infile_1 not required if --outdir\n") if($infile_1);
	die("--infile_2 not required if --outdir\n") if($infile_2);
	die("--prefix not required if --outdir\n") if($prefix);
}else{
	if(!$infile_1 || !$infile_2){
		die "--infile_1 and --infile_2 required.\n";
	}
	die("--infile_1 file is empty or does not exists! (Typed wrong filename?)\n") if((!-e $infile_1) and (!-s $infile_1));
	die("--infile_2 file is empty or does not exists! (Typed wrong filename?)\n") if((!-e $infile_2) and (!-s $infile_2));
}
die("--num_threads arg required\n") unless($num_threads);
die("no \$TMPDIR variable defined\n") unless($TMPDIR);
die("--outdir does not exists...\n") unless(-d $outdir);

my $tmpdir = File::Temp->newdir(
    "tmpDirFlashWrapperXXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 # 1-Delete after execution. 0-Do not delete.
);

$m = 10   unless($m);
$M = 280  unless($M);
$x = 0.25 unless($x);
$p = 33   unless($p);
$s = 5    unless($s);
$r = 150  unless($r);
$f = 250  unless($f);

## MAIN
$n = 20000 unless($n);
$start_at = 0 unless($start_at);
mkdir $outdir."/reads_subset" unless -d $outdir."/reads_subset";
mkdir $outdir."/assembly_complete" unless -d $outdir."/assembly_complete";
mkdir $outdir."/assembly_sampling" unless -d $outdir."/assembly_sampling";
mkdir $outdir."/Qscores" unless -d $outdir."/Qscores";

my %hash1;
my %hash2;
my %hashAssembled;
my %hashLog;

if($indir && $sampleSheet){
#	my $rHoAoH_sampleInfo  = SampleSheet::parseSampleSheetAsHash($sampleSheet);
#	
#	# Loop through each samples in sample sheet
#	my $rawReadDir    = $outdir."/raw_reads";
#	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
#		my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
#		my $R1_gz;
#		my $R2_gz;
#		my $R1_name;
#		my $R2_name;
#
#		for my $rH_laneInfo (@$rAoH_sampleLanes) { #this rH_laneInfo contains the complete line info from the sample sheet for this sample.
#			$R1_gz = $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};
#		 	$R2_gz = $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read2File'};
#			$R1_name = $rH_laneInfo->{'read1File'};
#		 	$R2_name = $rH_laneInfo->{'read2File'};
#		}
#
#		# Uncompress
#		$R1_name =~ s/\.gz//;
#		$R2_name =~ s/\.gz//;
#		my $R1_fastq = $tmpdir."/".$R1_name;
#		my $R2_fastq = $tmpdir."/".$R2_name;
#		system("gunzip -c ".$R1_gz." > ".$R1_fastq);
#		$? != 0 ? die "command failed: $!\n" : print STDERR $R1_gz." decompressed... into ".$R1_fastq."\n" if($verbose);
#		system("gunzip -c ".$R2_gz." > ".$R2_fastq);
#		$? != 0 ? die "command failed: $!\n" : print STDERR $R2_gz." decompressed... into ".$R2_fastq."\n" if($verbose);
#
#		print STDERR $R1_fastq."\n" if($verbose);
#		print STDERR $R2_fastq."\n" if($verbose);
#
#		# Assemble
#		my $assembledReads = assemble($R1_fastq, $R2_fastq, $outdir);
#		my $assembledName = basename($assembledReads);	
#		$assembledName =~ s/\.fastq//;
#
#		# Compute Qstats and store qstat sheet in a hash for later retrival when building qscore plots.
#		$hash1{'fastq'}{$R1_name} = $R1_fastq;
#		$hash2{'fastq'}{$R2_name} = $R2_fastq;
#		$hashAssembled{'fastq'}{$assembledName} = $assembledReads;
#
#	}

}elsif($infile_1 && $infile_2){
	assemble($infile_1, $infile_2, $outdir);
}

my ($qscores1, $qscores2, $qscoresAss) = generateQscoreSheet(\%hash1, \%hash2, \%hashAssembled) if($computeQscores);	
generateQscorePlots($qscores1, $qscores2, $qscoresAss) if($computeQscores);

writeLog();

exit;

## SUBROUTINES

sub generateQscorePlots{
	my $qs1 = shift;
	my $qs2 = shift;
	my $qsA = shift;

	# Generate plots for R1 and R2
	my $cmd = "~/build/tools/generateQscorePlots.pl";
	$cmd .= " --infile_1 ".$qs1;
	$cmd .= " --infile_2 ".$qs2;
	$cmd .= " --name R1_R2";
	$cmd .= " --pdf ".$outdir."/R1R2Qplots.pdf";
	$cmd .= " --display 1";
	$cmd .= " --paired";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDERR "R1R2 Qplots generated...\n" if($verbose);

	# Generate plots for assembled reads.
	$cmd = "~/build/tools/generateQscorePlots.pl";
	$cmd .= " --infile_1 ".$qsA;
	$cmd .= " --name Assembled_reads";
	$cmd .= " --pdf ".$outdir."/assembledReadsQplots.pdf";
	$cmd .= " --display 1";
	$cmd .= " --single";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDERR "Extended fragments Qplots generated...\n" if($verbose);
}

sub writeLog{

	open(OUT, '>'.$outdir."/summaryLog.txt") or die "Can't open file ".$outdir."/summaryLog.txt";
	print OUT "#SampleName\tTotal reads\tCombined reads\tUncombined reads\tPercent combined\n";

	for my $key ( keys %hashLog ) {
        my $file = $hashLog{$key};
		
		my $totalReads;
		my $combinedReads;
		my $uncombinedReads;
		my $percentCombined;
		
		open(IN, '<'.$file) or die "Can't open file ".$file."\n";
		while(<IN>){
			chomp;	
			if($_ =~ m/Total reads\:\s*(\d+)/){
				$totalReads = $1;
			}elsif($_ =~ m/Combined reads\:\s*(\d+)/){
				$combinedReads = $1;
			}elsif($_ =~ m/Uncombined reads\:\s*(\d+)/){
				$uncombinedReads = $1;
			}elsif($_ =~ m/Percent combined\:\s*(\d+\.\d+\%)/){
				$percentCombined = $1;
			}
		}
		print OUT "$key\t$totalReads\t$combinedReads\t$uncombinedReads\t$percentCombined\n";
    }
	close(OUT);
}

sub generateQscoreSheet{

	my( $href1, $href2, $hrefAssembled) = @_;
	my @prefixes =("R1", "R2", "Assembled");
	
	my @AoH = ($href1, $href2, $hrefAssembled);
	my @Qscores;	

	# Get file names (with full path).
	for my $href ( @AoH ) {
		my %hash = %$href;	
		my $prefix = shift @prefixes;
		my $files = "";
		foreach my $key (keys %hash) {
			foreach my $key2 (keys %{ $hash{$key} }) {
				$files .= " --fastq ".$hash{'fastq'}{$key2};
			}
			my $cmd = "~/build/tools/generateQscoreSheets.pl ";
			$cmd .= $files;
			$cmd .= " --outfile ".$outdir."/Qscores/qscore_".$prefix.".txt";
			$cmd .= " --phred 33";
			$cmd .= " --num_threads ".$num_threads;
			system($cmd);
			$? != 0 ? die "command failed: $!\n" : print STDERR "Reads subset generated...\n" if($verbose);
			push(@Qscores, $outdir."/Qscores/qscore_".$prefix.".txt");
	    }
	}
	return @Qscores;	
}

sub assemble{

	my $infile_1 = shift;
	my $infile_2 = shift;
	my $outdir = shift;

	my $R1Prefix = basename($infile_1);
	my $R2Prefix = basename($infile_2);
	$R1Prefix =~ s/\.pair1\.fastq//;	
	$R2Prefix =~ s/\.pair2\.fastq//;
	my $ASSPrefix = $R1Prefix;
	$ASSPrefix =~ s/\.pair1\.fastq//;
	$ASSPrefix =~ s/\_1\.fastq//;
	$ASSPrefix =~ s/\_R1\.fastq//;
	$ASSPrefix =~ s/\_read1\.fastq//;
	$ASSPrefix = $prefix if($prefix);

	# First get subset of reads (random).
	my $bin = `which getSubset.pl`; die if(!defined($bin)); chomp($bin);
	my $cmd = $bin;
	$cmd .= " --infile_1 ".$infile_1;
	$cmd .= " --infile_2 ".$infile_2;
	$cmd .= " --outfile_1 ".$outdir."/reads_subset/".$R1Prefix."_subset1.fastq";
	$cmd .= " --outfile_2 ".$outdir."/reads_subset/".$R2Prefix."_subset2.fastq";
	$cmd .= " --n ".$n;
	system($cmd) if($start_at <= 0);
	$? != 0 ? die "command failed: $!\n" : print STDERR "Reads subset generated...\n" if($verbose);
	
	# Then do assembly on subset and get stats for assembly.
	my ($averageReadLength, $averageFragLength, $stddev) = getStatsForAssembly($outdir."/reads_subset/".$R1Prefix."_subset1.fastq", $outdir."/reads_subset/".$R2Prefix."_subset2.fastq", $n);
	$averageReadLength = sprintf("%.0f", $averageReadLength);
	$averageFragLength = sprintf("%.0f", $averageFragLength);
	$stddev = sprintf("%.0f", $stddev);
	
	print STDERR "Average read length:\t".$averageReadLength."\nAverage fragment length:\t".$averageFragLength."\nStandard deviation:\t".$stddev."\n";	

	# Then do "real" assembly	
	$cmd = "flash ";
	$cmd .= " ".$infile_1;
	$cmd .= " ".$infile_2;
	$cmd .= " -r ".$averageReadLength;
	$cmd .= " -s ".$stddev;
	$cmd .= " -f ".$averageFragLength;
	$cmd .= " -o ".$ASSPrefix;
	$cmd .= " -d ".$outdir."/assembly_complete";
	$cmd .= " -t ".$num_threads;
	$cmd .= " &> ".$outdir."/assembly_complete/".$ASSPrefix.".log";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDERR "FLASH assembly for sampling done...\n" if($verbose);
	
	$hashLog{$ASSPrefix} = $outdir."/assembly_complete/".$ASSPrefix.".log";

	return $outdir."/assembly_complete/".$ASSPrefix.".extendedFrags.fastq";
}

# Remove temp files
sub END{
	local $?;
	system("rm ".$tmpdir." -rf");
}

# TODO Put in separate script.
# Return some paired-end overlapping assembly using subsampling of library
# $_[0] = reads1, $_[1] = reads2
# 1-read length, 2-fragment length, 3-stdeviation
sub getStatsForAssembly{
    my $reads_1         = shift;
    my $reads_2 	    = shift;
	my $sampling_number = shift;

    my @readLength = ();
    my @fragment_length = ();
    my $std = 0;

	# Reads 1
    my $counter = 0;
    my $in = new Iterator::FastqDb($reads_1) or die("Unable to open Fastq file, $reads_1\n");
    while( my $curr = $in->next_seq() and ($counter < $sampling_number) ){
        push(@readLength, length($curr->seq()));
        $counter++;
    }

	# Reads 2
    $counter = 0;
    $in = new Iterator::FastqDb($reads_2) or die("Unable to open Fastq file, $reads_2\n");
    while( my $curr = $in->next_seq() and ($counter < $sampling_number) ){
        push(@readLength, length($curr->seq()));
        $counter++;
    }
    my $averageReadLength = sum(@readLength)/@readLength;
    my $stat = Statistics::Descriptive::Sparse->new();

	# Then do assembly on subset and get stats for assembly.
	my $cmd = "flash ";
	$cmd .= " ".$reads_1;
	$cmd .= " ".$reads_2;
	$cmd .= " -r ".$averageReadLength;
	$cmd .= " -o assembled";
	$cmd .= " -d ".$outdir."/assembly_sampling";
	$cmd .= " -t ".$num_threads;
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDERR "FLASH assembly for sampling done...\n" if($verbose);

    $counter = 0;
    $in = new Iterator::FastqDb($outdir."/assembly_sampling/assembled.extendedFrags.fastq") or die("Unable to open Fastq file, ".$outdir."/assembly_sampling/assembled.extendedFrags.fastq");
    while( my $curr = $in->next_seq() and ($counter < $sampling_number) ){
        $stat->add_data(length($curr->seq()));
        push(@fragment_length, length($curr->seq()));
        $counter++;
    }
    my $averageFragLength = sum(@fragment_length)/@fragment_length;
    $std = $stat->standard_deviation();
    my $stdevOfFragLength = ($std/$averageFragLength)*100;
	$std = 1 if($std == 0); #Because FLASH does not handle stdev values of 0.	
    
	return($averageReadLength, $averageFragLength, $std);
}
