#!/usr/bin/env perl

use strict;
use warnings;

BEGIN{
	#Makesure we can find the GetConfig::LoadModules module relative to this script install
	use File::Basename;
	use Cwd 'abs_path';
	my ( undef, $mod_path, undef ) = fileparse( abs_path(__FILE__) );
	unshift @INC, $mod_path."lib";
	unshift @INC, $mod_path."../../lib";
}

## LIBRARIES
use Env qw/TMPDIR PATH PYTHONPATH PERL5LIB/;
use Getopt::Long;
use POSIX qw/strftime/;
use Sys::Hostname;
use Statistics::Descriptive;
use List::Util qw(sum);
use PDF::Create;
use File::Spec::Functions qw(rel2abs);
use File::Slurp;
use File::Temp;
use File::Which;
use File::Basename;

# Libraries mugqic pipeline
use RRNAAmplicons;
use MicrobialEcology;
use GqSeqUtils;
use LoadConfig;
use SampleSheet;
use SubmitToCluster;
use Job;
use Tools;
use Version;

our $VERSION = "0.5";
$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
RRNATagger.pl

PURPOSE:
Pipeline for rRNA 16S/18S/ITS amplicons.

Steps:
Steps will vary depending on the type of analysis being run (--lib_type <string> and --fungi_ITS, --bactArch, --everything).

INPUT:
--external_infile <string> : We assume it is a paired-end fastq library file. Can be .gz
OR
--sampleSheet <string>     : MUGQIC sample sheet.

--configFile <string>      : MUGQIC style configuration .ini file.
--barcodes <string>        : Barcode file in fasta format. Multiple barcode files can be supplied 
                             if for instance, this is a run having samples from multiple users.
--config_file <string>     : Configuration file containing paths to the RDP training set, contaminant
                             database, PhiX fasta reference and chimera fasta reference database.
                             This ini file also contains all the paramters for the various steps
                             of the pipeline. 

ONE of the three following flags:
--fungi_ITS                : Put this flag if analyzing fungi ITS.
--bactArch                 : If analyzing bacteria/archeal.
--everything               : If analyzing everything.

ADVANCED OPTIONS
--start_at <int>           : Pipeline will be start at step <int> (optional - mostly for debugging)
--end_at <int>             : Pipeline will stop at step <int> (optional - mostly for debugging)

OUTPUT:
--outdir <string>          : Output directory into which files generated through the pipeline will 
                             be stored.

NOTES:
This pipeline can also process pyrotag data, but it has to 
be submitted as properly formatted fastq files with barcodes 
in headers etc.

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
Also see Qiime tools - www.qiime.org	

ENDHERE

## OPTIONS
my (
	$help, $external_infile, $log, $prefix, @barcodes,
	$start_at, $end_at, $fungi_ITS, $bactArch, $everything, $no_clustering, $sampleSheet,
	$run_on_cluster, $config_file, $noMsub
);
my $verbose = 0;

my @cmd = @ARGV;

GetOptions(
  'external_infile=s' => \$external_infile,
	'sampleSheet=s'		  => \$sampleSheet,
	'barcodes=s' 			  => \@barcodes,
	'config_file=s'		  => \$config_file,
	'fungi_ITS'				  => \$fungi_ITS,
	'bactArch'          => \$bactArch,
	'everything'        => \$everything,
	'run_on_cluster' 	  => \$run_on_cluster,
	'start_at=i' 			  => \$start_at,
	'end_at=i' 				  => \$end_at,
	'no_clustering'     => \$no_clustering,
	'noMsub'				    => \$noMsub,
  'verbose' 			    => \$verbose,
  'help' 					    => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--fungi_ITS or --bactArch or --everything required\n") unless($fungi_ITS || $bactArch || $everything);
die("--config_file arg missing\n") unless($config_file);
my $rHoAoH_sampleInfo;
if($external_infile){
    die("--external_infile does not exists! Typed wrong filename?\n") if((!-e $external_infile) and (!-s $external_infile));
	die "Please do not supply a sample sheet if you are supplying an external fastq library...\n" if($sampleSheet) ;
}else{
	$rHoAoH_sampleInfo  = SampleSheet::parseSampleSheetAsHash($sampleSheet);
}

die("--config_file does not exists! Typed wrong filename?\n") if((!-e $config_file) and (!-s $config_file));

my %cfg               = LoadConfig->readConfigFile($config_file);
my $iniFileTmpDir     = LoadConfig::getParam(\%cfg, 'default', 'tmpDir');
my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerXXXXXXX",
    DIR => $iniFileTmpDir."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

mkdir $tmpdir unless -d $tmpdir;

my $outdir             = LoadConfig::getParam(\%cfg, 'default', 'currentDir');
my $clusteringMethod   = LoadConfig::getParam(\%cfg, 'clustering', 'clusteringMethod'); 
my $projectName        = LoadConfig::getParam(\%cfg, 'default', 'projectName');
my $lib_type           = LoadConfig::getParam(\%cfg, 'default', 'libraryType');
$projectName           = "RRNATagger-project" unless($projectName); 

print STDERR "[DEBUG] **************************************************************** \n";
print STDERR "[DEBUG] ********* Generating RRNATagger pipeline commands... *********** \n";
print STDERR "[DEBUG] **************************************************************** \n";

die("--outdir required\n") unless $outdir;
die("--barcodes file containing barcodes in fasta format required\n") unless @barcodes > 0;
print STDERR "[DEBUG] Using clustering method ".$clusteringMethod."\n";	
foreach my $barcodes (@barcodes){
	print STDERR "[DEBUG] Using barcodes: ".$barcodes."\n" if($verbose);
	die("--barcodes file is empty or does not exists! (Typed wrong filename?)\n") if((!-e $barcodes) and (!-s $barcodes));
}

## Place current dir in $ENV{PATH} values
my $bin_dir_name = dirname(rel2abs($0));
$ENV{PATH} = $bin_dir_name.":".$ENV{PATH};

die("--lib_type: Please enter 'ex' or ",
	"'nc1nc2' or",
	"'nc1'\n") if(
		$lib_type ne "ex",
		and $lib_type ne "nc1nc2",
		and $lib_type ne "nc1",
);

my $fwd_primer  = LoadConfig::getParam(\%cfg, 'itags_QC', 'forwardPrimer', 0, 'filepath');
my $rev_primer  = LoadConfig::getParam(\%cfg, 'itags_QC', 'reversePrimer', 0, 'filepath');

$start_at = 1 unless $start_at;
$end_at = 999 unless $end_at;
my $currStep = 1;
$outdir = substr($outdir, 0, length($outdir)-1 ) if substr($outdir, -1) eq "/";
mkdir $outdir unless -d $outdir;
my $compute_cluster_outdir = $outdir."/compute_cluster_outfiles";
mkdir $compute_cluster_outdir unless -d $compute_cluster_outdir;

SubmitToCluster::initPipeline;

#print STDOUT "export PATH=./scripts:\$PATH\n";

my $infile;
my $dependency = undef;
if($external_infile){
	$infile = $external_infile; # If reads are coming from elsewhere than GQIC
}elsif($sampleSheet){
	
	if($start_at <= $currStep && $end_at >= $currStep){
		# Concatenate all libraries. Loop through .fastq.gz that will be specified by sampleSetup.pl.
		# Basically, concatenate all *R1.fastq.gz and *R2.fastq.gz separately. Assuming barcodes
		# are in headers.
		# Loop through each samples in sample sheet
		my $rO_jobMergeBarcodes = RRNAAmplicons::mergeBarcodes(
			\%cfg,
			$rHoAoH_sampleInfo,
			$outdir
		);
		if(!$rO_jobMergeBarcodes->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "mergeBarcodes", "mergeBarcodes" , 'MERGEBARCODES', $dependency,  "global", $rO_jobMergeBarcodes) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobMergeBarcodes->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
	}
	
	$infile = $outdir."/raw_reads/paired.fastq";
}else{
	die "--sampleSheet or --external_infile arg missing...\n"
}
print STDERR "[DEBUG] Step #".$currStep." mergeBarcodes\n";
$currStep++;

## MAKE COMMON OUTDIR FOR A LIBRARY
my $duk_dir = $outdir."/duk";
print STDOUT "mkdir -p ".$duk_dir."\n";
system("mkdir -p ".$duk_dir);

## RUN PIPELINE FOR EACH BARCODES FILES
my $loop_counter = 1;
my $curr_dir;
my @sr_array; # For single reads pipeline

# number of samples will be stored in this following variable:
my $numberOfSamples = 0;

foreach my $barcodes (@barcodes){

  ## Count samples being processed.
  my $currNumberOfSamples = `grep -o ">" $barcodes | wc -l`;
  chomp($currNumberOfSamples);
  $numberOfSamples += $currNumberOfSamples;
  print STDERR "[DEBUG] Number of samples: $numberOfSamples\n";

	# Declare variables holding file name of reads going into clustering.
	# These will be QC passed reads.	
	my $assembled_filtered;
	my $reads_1_filtered;
	my $reads_2_filtered;
	
	# Make barcodes specific dirs
	my $barcodes_filename = basename($barcodes);
	
	# Try to remove fasta extensions.
	if($barcodes_filename =~ m/(.*)\.fasta|fsa|fa|fas|fna/i){
		$barcodes_filename = $1;
	}	
	$curr_dir = $outdir."/".$barcodes_filename;	
	print STDOUT "mkdir -p ".$curr_dir."\n";	
	system("mkdir -p ".$curr_dir);	
	
	# In each barcodes specific dirs, make a fastqs dir in which will be the trimmed filtered, assembled, etc. fastqs...	
	print STDOUT "mkdir -p ".$curr_dir."/fastqs\n";
	system("mkdir -p ".$curr_dir."/fastqs");

	# Maybe eventually change that to give more meaningful names...
	my $assembled_dir = $outdir."/".$barcodes_filename."/assembled";
	my $reads_1_dir   = $outdir."/".$barcodes_filename."/reads_1";
	my $reads_2_dir   = $outdir."/".$barcodes_filename."/reads_2";
	
	# MAKE DIRS
	make_dirs($assembled_dir) if($lib_type eq "exnc1nc2" ,
								or $lib_type eq "ex" 
							);
	make_dirs($reads_1_dir) if($lib_type eq "exnc1nc2",
								or $lib_type eq "nc1"
							);
	make_dirs($reads_2_dir) if($lib_type eq "nc2",
								or $lib_type eq "nc1nc2",
								or $lib_type eq "exnc1nc2"
							);

	## MAIN 

	# START EXECUTION
	my $infile_gzip; 
	my $infile_fastq;
	if($infile =~ m/\.gz|\.gzip/){
		$infile_gzip = $infile;
	}elsif($infile =~ m/\.fq|\.fastq/){
		$infile_fastq = $infile;
	}
	
	# Remove contaminants onluy once for ressources savings.
	if($loop_counter == 1){
		# Only perform contaminants filtering one time
		# Remove contam reads
		if($infile_fastq){
			my $rO_jobDukWrapper = RRNAAmplicons::dukWrapper(
				\%cfg,
				$infile_fastq,
				$duk_dir."/contam.fastq",
				$duk_dir."/ncontam.fastq",
				$duk_dir."/duk_contam_log.txt",
				LoadConfig::getParam(\%cfg, 'DB',  'contaminants')
			);
			if(!$rO_jobDukWrapper->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "duk_wrapper", "duk_wrapper_contam"."_".$loop_counter, 'DUKWRAPPERCONTAM'."_".$loop_counter, $dependency,  "global", $rO_jobDukWrapper) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobDukWrapper->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\tStep #".$currStep." duk_wrapper_contam\n";
			$currStep++;

		}elsif($infile_gzip){
			my $rO_jobDuk = RRNAAmplicons::duk(
				\%cfg,
				$duk_dir."/duk_contam_log.txt",
				$duk_dir."/ncontam.fastq",
				$duk_dir."/contam.fastq",
				LoadConfig::getParam(\%cfg, 'DB',  'contaminants'),
				$infile_gzip
			);
			if(!$rO_jobDuk->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "duk", "duk_contam"."_".$loop_counter , 'DUKCONTAM'."_".$loop_counter, $dependency,  "global", $rO_jobDuk) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobDuk->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\tStep #".$currStep." duk_contam\n";
			$currStep++;
		}
	}
	
	# write to log for each barcodes but just do this step once for cpu usage savings.
	if($loop_counter == 1){
		my $rO_jobDukWrapperPhix = RRNAAmplicons::dukWrapper(
			\%cfg,
			$duk_dir."/ncontam.fastq",
			$duk_dir."/ncontam_phix.fastq",
			$duk_dir."/ncontam_nphix.fastq",
			$duk_dir."/duk_phix_log.txt",
			LoadConfig::getParam(\%cfg, 'DB',  'contaminants')
		);
		if(!$rO_jobDukWrapperPhix->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "duk_wrapper", "duk_wrapper_phix"."_".$loop_counter , 'DUKWRAPPERPHIX'."_".$loop_counter, $dependency,  "global", $rO_jobDukWrapperPhix) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobDukWrapperPhix->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." duk_wrapper_phix\n";
		$currStep++;
	}
	
	# Before moving foward, split the library upon relevant barcodes. one file having all good barcodes
	# important in case one MiSeq/HiSeq lane contains samples from more than one user/project.
	my $barcodesOutput;
	if($lib_type eq "nc1"){
		$barcodesOutput = $curr_dir."/fastqs/ncontam_nphix_1.fastq"	
	}else{
		$barcodesOutput = $curr_dir."/fastqs/ncontam_nphix.fastq"		
	}
	my $rO_jobSplitBarcodes = RRNAAmplicons::splitBarcodes(
		\%cfg,
		$duk_dir."/ncontam_nphix.fastq",
		$barcodes,
		$barcodesOutput,
		$curr_dir."/fastqs/ncontam_nphix_barcodes.txt"
	);
	if(!$rO_jobSplitBarcodes->isUp2Date()) {
		SubmitToCluster::printSubmitCmd(\%cfg, "barcodes", "barcodes"."_".$loop_counter , 'SPLITBARCODES'."_".$loop_counter, $dependency,  "global", $rO_jobSplitBarcodes) if($start_at <= $currStep && $end_at >= $currStep); 
		$dependency = $rO_jobSplitBarcodes->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
	}
	print STDERR "[DEBUG]\tStep #".$currStep." barcodes\n";
	$currStep++;

	if($lib_type eq "ex"){
		# Remove unpaired reads
		my $rO_jobRemoveUnpairedReads = RRNAAmplicons::removeUnpairedReads(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix.fastq",
			$curr_dir."/fastqs/ncontam_nphix_paired.fastq",
			$curr_dir."/fastqs/ncontam_nphix_unpaired_1.fastq",
			$curr_dir."/fastqs/ncontam_nphix_unpaired_2.fastq"
		);
		if(!$rO_jobRemoveUnpairedReads->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "remove_unpaired", "remove_unpaired"."_".$loop_counter , 'REMOVEUNPAIRED'."_".$loop_counter, $dependency, "global", $rO_jobRemoveUnpairedReads) if($start_at <= $currStep && $end_at >= $currStep); 
		}
		$dependency = $rO_jobRemoveUnpairedReads->getCommandJobId(0);
		print STDERR "[DEBUG]\tStep #".$currStep." remove_unpaired\n";
		$currStep++;
		
		# Split paired reads in reads_1 and reads_2.
		my $rO_jobSplitPairs = RRNAAmplicons::splitPairs(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix_paired.fastq",
			$curr_dir."/fastqs/ncontam_nphix_1.fastq",
			$curr_dir."/fastqs/ncontam_nphix_2.fastq"
		);
		if(!$rO_jobSplitPairs->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "split_pairs", "split_pairs"."_".$loop_counter, "SPLITPAIRS"."_".$loop_counter, $dependency, "global", $rO_jobSplitPairs) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobSplitPairs->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." split_pairs\n";
		$currStep++;
	}

		
	if($lib_type =~ m/nc1/ || $lib_type =~ m/ex/){
		# generate qscores and store them in only one sheet. Reads 1
		my $rO_jobQscoreSheetsR1 = RRNAAmplicons::generateQscoreSheet(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix_1.fastq",
			"reads_1",
			$curr_dir."/fastqs/qual_stats_1_log.txt",
			$curr_dir."/fastqs/qscores_1.tsv",
			$barcodes
		);
		if(!$rO_jobQscoreSheetsR1->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_sheet", "qscore_sheet_R1_raw"."_".$loop_counter, "QSCORESHEETR1"."_".$loop_counter, $dependency, "global", $rO_jobQscoreSheetsR1) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreSheetsR1->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_sheet_R1_raw\n";
		$currStep++;
	}	

	if($lib_type =~ m/nc2/ || $lib_type =~ m/ex/){
		# generate qscores and store them in only one sheet. Reads 2
		my $rO_jobQscoreSheetsR2 = RRNAAmplicons::generateQscoreSheet(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix_2.fastq",
			"reads_2",
			$curr_dir."/fastqs/qual_stats_2_log.txt",
			$curr_dir."/fastqs/qscores_2.tsv",
			$barcodes
		);
		if(!$rO_jobQscoreSheetsR2->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_sheet", "qscore_sheet_R2_raw"."_".$loop_counter, "QSCORESHEETR2"."_".$loop_counter, $dependency, "global", $rO_jobQscoreSheetsR2) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreSheetsR2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_sheet_R2_raw\n";
		$currStep++;
	}

	if($lib_type =~ m/ex/ || ($lib_type =~ m/nc1/ && $lib_type =~ m/nc2/) ){
		# generate pdf graph.
		my $rO_jobQscoreGraphR1R2 = RRNAAmplicons::generateQscoreGraphPaired(
			\%cfg,
			$curr_dir."/fastqs/qscores_1.tsv",
			$curr_dir."/fastqs/qscores_2.tsv",
			$curr_dir."/qual_stats_unassembled.pdf"
		);
		if(!$rO_jobQscoreGraphR1R2->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_plot", "qscore_plot_R1R2_raw"."_".$loop_counter, "QSCOREPLOTR1R2RAW"."_".$loop_counter, $dependency, "global", $rO_jobQscoreGraphR1R2) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreGraphR1R2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_plot_R1R2_raw\n";
		$currStep++;

	}elsif($lib_type eq "nc1"){
		# generate pdf graph.
		my $rO_jobQscoreGraphR1 = RRNAAmplicons::generateQscoreGraphSingle(
			\%cfg,
			$curr_dir."/fastqs/qscores_1.tsv",
			"qual_stats_unassembled",
			$curr_dir."/qual_stats_unassembled.pdf"
		);
		if(!$rO_jobQscoreGraphR1->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_plot", "qscore_plot_R1_raw"."_".$loop_counter, "QSCOREPLOTR1RAW"."_".$loop_counter, $dependency, "global", $rO_jobQscoreGraphR1) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreGraphR1->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_plot_R1_raw\n";
		$currStep++;

	}elsif($lib_type eq "nc2"){
		# generate pdf graph.
		my $rO_jobQscoreGraphR2 = RRNAAmplicons::generateQscoreGraphSingle(
			\%cfg,
			$curr_dir."/fastqs/qscores_2.tsv",
			"qual_stats_unassembled",
			$curr_dir."/qual_stats_unassembled.pdf"
		);
		if(!$rO_jobQscoreGraphR2->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_plot", "qscore_plot_R2_raw"."_".$loop_counter, "QSCOREPLOTR2RAW"."_".$loop_counter, $dependency, "global", $rO_jobQscoreGraphR2) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreGraphR2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_plot_R2\n";
		$currStep++;
	}	

	# Cut reads prior to assembly.	
	if($lib_type eq "ex"){	
		
	    # This step is essentially to prepare reads for paired-end overlapping assembly. I.e. 2x250 reads needs to be shorten to about 165 nt to get proper assembly.
	    # If 2x150, the first and last one or two nt are often of bad qual. So just trim them. Skip this step only for reads_1_only_paired_library.	
		
		# Cut reads 1
		my $rO_jobCutR1 = RRNAAmplicons::cutReads(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix_1.fastq",
			LoadConfig::getParam(\%cfg, 'itags_QC', 'R1_start'),
			(LoadConfig::getParam(\%cfg, 'default',  'readLength') - LoadConfig::getParam(\%cfg, 'itags_QC', 'R1_end')),
			$curr_dir."/fastqs/ncontam_nphix_trimmed_1.fastq"
		);
		if(!$rO_jobCutR1->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "trim", "trim_R1"."_".$loop_counter, "TRIMR1"."_".$loop_counter, $dependency, "global", $rO_jobCutR1) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobCutR1->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." trim_R1\n";
		$currStep++;

		# Cut reads 2
		my $rO_jobCutR2 = RRNAAmplicons::cutReads(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix_2.fastq",
			LoadConfig::getParam(\%cfg, 'itags_QC', 'R2_start'),
			(LoadConfig::getParam(\%cfg, 'default',  'readLength') - LoadConfig::getParam(\%cfg, 'itags_QC', 'R2_end')),
			$curr_dir."/fastqs/ncontam_nphix_trimmed_2.fastq"
		);
		if(!$rO_jobCutR2->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "trim", "trim_R2"."_".$loop_counter, "TRIMR2"."_".$loop_counter, $dependency, "global", $rO_jobCutR2) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobCutR2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." trim_R2\n";
		$currStep++;
		
		# flash wrapper which will Subsample lib for flash
		my $rO_jobFlash = RRNAAmplicons::flash(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix_trimmed_1.fastq",
			$curr_dir."/fastqs/ncontam_nphix_trimmed_2.fastq",
			"ncontam_nphix_trimmed",
			$curr_dir."/fastqs/"
		);
		if(!$rO_jobFlash->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "flash", "flash"."_".$loop_counter, "FLASH"."_".$loop_counter, $dependency, "global", $rO_jobFlash) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobFlash->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." flash\n";
		$currStep++;

		# Remove primers (if primers appear in reads) and do Itags QC (filter reads).
		#if($rev_primer || $fwd_primer ){ #Don't execute this step if primers are not to be removed.	
			#if(defined $rev_primer && defined $fwd_primer){
      #  # Leave as is...
			#}elsif(!defined $fwd_primer && !defined $rev_primer){
			#	$fwd_primer = "null";	
			#	$rev_primer = "null";	
			#}else{
      #  if(!defined $fwd_primer){
			#	  $fwd_primer = "null";
      #  }	
      #  if(!defined $rev_primer){
			#	  $rev_primer = "null";
      #  }	
			#}
      $fwd_primer = "null" if(!defined $fwd_primer);
      $rev_primer = "null" if(!defined $rev_primer);
      $fwd_primer = "null" if($fwd_primer eq "");
      $rev_primer = "null" if($rev_primer eq "");
      print STDERR "[DEBUG] fwd_primer: ".$fwd_primer."\n";	
      print STDERR "[DEBUG] fwd_primer: ".$rev_primer."\n";	
			my $rO_jobItagsQC = RRNAAmplicons::itagsQC(
				\%cfg,
				$curr_dir."/fastqs/assembly_complete/ncontam_nphix_trimmed.extendedFrags.fastq",
				$rev_primer,
				$fwd_primer,
				$curr_dir."/fastqs/ncontam_nphix_trimmed.extendedFrags_QCpassed.fastq",
				$curr_dir."/fastqs/ncontam_nphix_trimmed.extendedFrags_QCfailed.fastq"	
			);
			if(!$rO_jobItagsQC->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "itags_QC", "itags_QC"."_".$loop_counter, "ITAGSQC"."_".$loop_counter, $dependency, "global", $rO_jobItagsQC) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobItagsQC->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\tStep #".$currStep." itags_QC\n";
			$currStep++;
		#}
		
		my $QC_passed_reads = $curr_dir."/fastqs/ncontam_nphix_trimmed.extendedFrags_QCpassed.fastq";
		$assembled_filtered = $QC_passed_reads;
			
		# Qual score graphs assembled + filtered reads
		my $rO_jobQscoreSheetsAss = RRNAAmplicons::generateQscoreSheet(
			\%cfg,
			$QC_passed_reads,
			"assembled_filtered",
			$curr_dir."/assembled/qscores/barcodes_assembled_QCpassed_log.txt",
			$curr_dir."/assembled/qscores/qscores_assembled_QCpassed.tsv",
			$barcodes
		);
		if(!$rO_jobQscoreSheetsAss->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_sheet", "qscore_sheet_ASS"."_".$loop_counter, "QSCORESHEETASS"."_".$loop_counter, $dependency, "global", $rO_jobQscoreSheetsAss) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreSheetsAss->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_sheet_ASS\n";
		$currStep++;
	
		# Generate pdf graph for assembled+filtered reads.
		my $rO_jobQscoreGraphAss = RRNAAmplicons::generateQscoreGraphSingle(
			\%cfg,
			$curr_dir."/assembled/qscores/qscores_assembled_QCpassed.tsv",
			"qual_stats_assembled_filtered",
			$curr_dir."/assembled/qscores/qscores_assembled_QCpassed.pdf"
		);
		if(!$rO_jobQscoreGraphAss->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_plot", "qscore_plot_ASS"."_".$loop_counter, "QSCOREPLOTASS"."_".$loop_counter, $dependency, "global", $rO_jobQscoreGraphAss) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreGraphAss->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_plot_ASS\n";
		$currStep++;
	}

	# Do reads 1
	if($lib_type =~ m/nc1/){
		
		# Remove primers (if primers appears in reads)
		my $currFwdPrimer;
		my $currRevPrimer;
		if($fwd_primer){
			$currFwdPrimer = $fwd_primer;
			$currRevPrimer = "null";	
		}else{
			$currFwdPrimer = "null";
			$currRevPrimer = "null";	
		}
		my $rO_jobItagsQCR1 = RRNAAmplicons::itagsQC(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix_1.fastq",
			$currRevPrimer,
			$currFwdPrimer,
			$curr_dir."/fastqs/ncontam_nphix_1_QCpassed.fastq",
			$curr_dir."/fastqs/ncontam_nphix_1_QCfailed.fastq"
		);
		if(!$rO_jobItagsQCR1->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "itags_QC", "itagsQCR1", "ITAGSQCR1"."_".$loop_counter, $dependency, "global", $rO_jobItagsQCR1) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobItagsQCR1->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." itags_QC\n";
		$currStep++;

		# Depending on what QC was done, put proper variables values.
		
		my $QC_passed_reads = $curr_dir."/fastqs/ncontam_nphix_1_QCpassed.fastq";
		$reads_1_filtered = $QC_passed_reads;
	
		# Qual score graphs for filtered reads
		my $rO_jobQscoreSheetsR1 = RRNAAmplicons::generateQscoreSheet(
			\%cfg,
			$QC_passed_reads,
			"reads_1_QCed",
			$curr_dir."/fastqs/barcodes_log_reads_1_QCpassed.txt",
			$curr_dir."/fastqs/qscores_reads_1_QCpassed.tsv",
			$barcodes
		);
		if(!$rO_jobQscoreSheetsR1->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_sheet", "qscore_sheet_R1QCed", "QSCORESHEETR1QCED"."_".$loop_counter, $dependency, "global", $rO_jobQscoreSheetsR1) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreSheetsR1->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_sheet_R1QCed\n";
		$currStep++;
	
		# Generate pdf graph for filtered reads 1
		my $rO_jobQscoreGraphR1QCed = RRNAAmplicons::generateQscoreGraphSingle(
			\%cfg,
			$curr_dir."/fastqs/qscores_reads_1_QCpassed.tsv",
			"qual_stats_reads_1_QCed",
			$curr_dir."/qscores_reads_1_QCpassed.pdf"
		);
		if(!$rO_jobQscoreGraphR1QCed->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_plot", "qscore_plot_R1QCed", "QSCOREPLOTR1QCED"."_".$loop_counter, $dependency, "global", $rO_jobQscoreGraphR1QCed) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreGraphR1QCed->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_plot_R1QCed\n";
		$currStep++;
	}
		
	# Do reads 2
	if($lib_type =~ m/nc2/){
		# Remove primers (if primers appears in reads)
		my $currFwdPrimer;
		my $currRevPrimer;
		if(defined $rev_primer){
			$currFwdPrimer = "null";
			$currRevPrimer = $rev_primer;	
		}else{
			$currFwdPrimer = "null";
			$currRevPrimer = "null";	
		}
		my $rO_jobItagsQCR2 = RRNAAmplicons::itagsQC(
			\%cfg,
			$curr_dir."/fastqs/ncontam_nphix_2.fastq",
			$currRevPrimer,
			$currFwdPrimer,
			$curr_dir."/fastqs/ncontam_nphix_2_QCpassed.fastq",
			$curr_dir."/fastqs/ncontam_nphix_2_QCfailed.fastq"
		);
		if(!$rO_jobItagsQCR2->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "itags_QC", "itagsQCR2", "ITAGSQCR2"."_".$loop_counter, $dependency, "global", $rO_jobItagsQCR2) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobItagsQCR2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." itagsQCR2\n";
		$currStep++;

		my $QC_passed_reads = $curr_dir."/fastqs/ncontam_nphix_2_noprimer_QCpassed.fastq";
		$reads_2_filtered = $QC_passed_reads;
	
		# Qual score graphs for filtered reads 1
		my $rO_jobQscoreSheetsR2 = RRNAAmplicons::generateQscoreSheet(
			\%cfg,
			$QC_passed_reads,
			"reads_2_QCed",
			$curr_dir."/fastqs/barcodes_log_reads_2_QCpassed.txt",
			$curr_dir."/fastqs/qscores_reads_2_QCpassed.tsv",
			$barcodes
		);
		if(!$rO_jobQscoreSheetsR2->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_sheet", "qscore_sheet_R2QCed", "QSCORESHEETR2QCED"."_".$loop_counter, $dependency, "global", $rO_jobQscoreSheetsR2) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreSheetsR2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_sheet_R2QCed\n";
		$currStep++;
	
		# Generate pdf graph for filtered reads 2
		my $rO_jobQscoreGraphR2QCed = RRNAAmplicons::generateQscoreGraphSingle(
			\%cfg,
			$curr_dir."/fastqs/qscores_reads_2_QCpassed.tsv",
			"qual_stats_reads_2_QCed",
			$curr_dir."/qscores_reads_2_QCpassed.pdf"
		);
		if(!$rO_jobQscoreGraphR2QCed->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "qscore_plot", "qscore_plot_R2QCed", "QSCOREPLOTR2QCED"."_".$loop_counter, $dependency, "global", $rO_jobQscoreGraphR2QCed) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobQscoreGraphR2QCed->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." qscore_plot_R2QCed\n";
		$currStep++;
	}

	# Prepare what reads will be clustered

	# Workflows (Reads clustering + classification) will be performed according to the value of $lib_type.
	if($lib_type =~ m/ex/ ){
		$dependency = executeWorkflow($assembled_dir, $assembled_filtered, $barcodes, "_overlapping", $dependency);
	}
	if($lib_type =~ m/nc1/){
		$dependency = executeWorkflow($reads_1_dir, $reads_1_filtered, $barcodes, "_reads1", $dependency);
	}
	if($lib_type =~  m/nc2/){
		$dependency = executeWorkflow($reads_1_dir, $reads_1_filtered, $barcodes, "_reads1", $dependency);
	}
	
	# Print count report
	my @reportFiles;
	my @reportNames;
	my $analysisType;
	my $barcodesDist;
	my $OTUtable;
	my $obsTable;
	
	push(@reportFiles, $infile);
	push(@reportFiles, $duk_dir."/contam.fastq");
	push(@reportFiles, $duk_dir."/ncontam_phix.fastq");
	push(@reportFiles, $duk_dir."/ncontam_nphix.fastq");
	push(@reportFiles, $curr_dir."/fastqs/ncontam_nphix_1.fastq");
	push(@reportFiles, $curr_dir."/fastqs/ncontam_nphix_2.fastq");
	push(@reportNames, "total_reads");
	push(@reportNames, "contaminants_reads");
	push(@reportNames, "phix_reads");
	push(@reportNames, "non_contam_non_phix_reads");
	push(@reportNames, "non_contam_non_phix_reads_1");
	push(@reportNames, "non_contam_non_phix_reads_2");
	
	if($lib_type =~ m/ex/){	
		push(@reportFiles, $curr_dir."/fastqs/assembly_complete/ncontam_nphix_trimmed.extendedFrags.fastq");	
		push(@reportFiles, $assembled_filtered);	
		push(@reportNames, "assembled_reads");	
		push(@reportNames, "assembled_reads_QC_passed");	
	}
	if($lib_type =~ m/nc1/){	
		push(@reportFiles, $reads_1_filtered);	
		push(@reportNames, "reads_1_QC_passed");
	}
	if($lib_type =~ m/nc2/){	
		push(@reportFiles, $reads_2_filtered );	
		push(@reportNames, "reads_2_QC_passed");
	}
	if($lib_type =~ m/ex/){	
		$analysisType = "assembled_clustered";
		$barcodesDist = $curr_dir."/fastqs/ncontam_nphix_barcodes.txt";;
		$OTUtable = $curr_dir."/assembled/otu_tables/otu_table.tsv";
		$obsTable = $curr_dir."/assembled/obs/obs_filtered.tsv";
	}
	if($lib_type =~ m/nc1/){	
		$analysisType = "reads_1_clustered";
		$barcodesDist = $curr_dir."/fastqs/barcodes_log_reads_1_QCpassed.txt";
		$OTUtable = $curr_dir."/reads_1/otu_tables/otu_table.tsv";
		$obsTable = $curr_dir."/reads_1/obs/obs_filtered.tsv";
	}
	if($lib_type =~ m/nc2/){	
		$analysisType = "reads_2_clustered";
		$barcodesDist = $curr_dir."/fastqs/barcodes_log_reads_2_QCpassed.txt";
		$OTUtable = $curr_dir."/reads_2/otu_tables/otu_table.tsv";
		$obsTable = $curr_dir."/reads_2/obs/obs_filtered.tsv";
	}

	my $rO_jobCountReport = RRNAAmplicons::countReport(
		\%cfg,
		\@reportFiles,
		\@reportNames,
		$analysisType,
		$barcodesDist,
		$OTUtable,
		$obsTable,
		$curr_dir."/countReport.tsv"
	);
	if(!$rO_jobCountReport->isUp2Date()) {
		SubmitToCluster::printSubmitCmd(\%cfg, "count_report", "count_report", "COUNTREPORT"."_".$loop_counter, $dependency, "global", $rO_jobCountReport) if($start_at <= $currStep && $end_at >= $currStep); 
		$dependency = $rO_jobCountReport->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
	}
	print STDERR "[DEBUG]\tStep #".$currStep." count_report\n";
	$currStep++;
	
	# Convert txt to pdf
	my $rO_jobTxtToPdf = RRNAAmplicons::txtToPdf(
		\%cfg,
		$curr_dir."/countReport.tsv",
		$curr_dir."/countReport.pdf"
	);
	if(!$rO_jobTxtToPdf->isUp2Date()) {
		SubmitToCluster::printSubmitCmd(\%cfg, "txtToPdf", "txtToPdf", "TXTTOPDF"."_".$loop_counter, $dependency, "global", $rO_jobTxtToPdf) if($start_at <= $currStep && $end_at >= $currStep); 
		$dependency = $rO_jobTxtToPdf->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
	}
	print STDERR "[DEBUG]\tStep #".$currStep." txtToPdf\n";
	$currStep++;
	
	my $qscore_pdf = $curr_dir."/qual_stats_unassembled.pdf";

	# R1, R2 and Assembled
	if($lib_type eq "exnc1nc2"){	
		mergePdf(
			$dependency,
			"exnc1nc2",
			$curr_dir."/log_barcodes_".$barcodes_filename."_count.pdf",
			$curr_dir."/reads_1/tax_summary/taxonomy_phylum_L2.pdf",
			$curr_dir."/reads_2/tax_summary/taxonomy_phylum_L2.pdf",
			$curr_dir."/assembled/tax_summary/taxonomy_phylum_L2.pdf",
			$qscore_pdf,
			$curr_dir."/REPORT_".$barcodes_filename.".pdf"
		);

	# R1 and R2
	}elsif($lib_type eq "nc1nc2"){
		mergePdf(
			$dependency,
			"nc1nc2",
			$curr_dir."/log_barcodes_".$barcodes_filename."_count.pdf",
			$curr_dir."/reads_1/tax_summary/taxonomy_phylum_L2.pdf",
			$curr_dir."/reads_2/tax_summary/taxonomy_phylum_L2.pdf",
			$qscore_pdf,
			$curr_dir."/REPORT_".$barcodes_filename.".pdf"
		);
	# R1 only
	}elsif($lib_type eq "nc1"){
		mergePdf(
			$dependency,
			"nc1",
			#$curr_dir."/log_barcodes_".$barcodes_filename."_count.pdf",
			$curr_dir."/countReport.pdf",
			$curr_dir."/reads_1/tax_summary/taxonomy_phylum_L2.pdf",
			$curr_dir."/reads_1/tax_summary/taxonomy_phylum_L2_relative.pdf",
			$qscore_pdf,
			$curr_dir."/REPORT_".$barcodes_filename.".pdf"
		);
	# or else...
		#my $prefix = shift(@_);
		#my $file_1 = shift(@_);
		#my $outfile = pop(@_);
	}elsif($lib_type eq "ex"){
		mergePdf(
			$dependency,
			"ex",
			$curr_dir."/countReport.pdf",
			$curr_dir."/assembled/tax_summary/taxonomy_phylum_L2.pdf",
			$curr_dir."/assembled/tax_summary/taxonomy_phylum_L2_relative.pdf",
			#$curr_dir."/assembled/tax_summary/taxonomy_phylum_bacteria_archaea_L2.pdf",
			$qscore_pdf,
			$curr_dir."/REPORT_".$barcodes_filename.".pdf"
		);
	}else{
		die "no valid \$lib_type value...\n";
	}
	
	#####################################################
	# Report / Deliverables
	#

	# Generate report with Noozle
	my $date = strftime "%Y-%m-%d", localtime;
	my $rO_jobDeliverables = RRNAAmplicons::clientReport(
		\%cfg,
		$config_file,
		$curr_dir,
		"RRNATagger",
		$curr_dir."/".$date
	);
	if(!$rO_jobDeliverables->isUp2Date()) {
		SubmitToCluster::printSubmitCmd(\%cfg, "report", "report", "NOZZLEREPORT"."_".$loop_counter, $dependency, "global", $rO_jobDeliverables) if($start_at <= $currStep && $end_at >= $currStep); 
		$dependency = $rO_jobDeliverables->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
	}
	print STDERR "[DEBUG]\tStep #".$currStep." Nozzle deliverables\n";
	$currStep++;
	
	$loop_counter++;
	## END COMMANDS
	
	# INNER SUBROUTINES START

	# Workflow in which a fastq library is clustered. Followed by diversity metric computation and taxonomic summary
	# @inputs = directory and fastq infile
	# @returns
	sub executeWorkflow{
		my($dir, $fastq, $barcodes, $prefix, $dependency) = @_;	
		
		if($clusteringMethod == 1){
			my $rO_jobClustering1 = RRNAAmplicons::clustering1(
				\%cfg,
				$fastq,
				$barcodes,
				$dir."/obs"
			);
			if(!$rO_jobClustering1->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "clustering", "clustering1", "CLUSTERING1"."_".$loop_counter, $dependency, "global", $rO_jobClustering1) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobClustering1->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." clustering_method_1\n";
			$currStep++;
		
		}elsif($clusteringMethod == 2){
			my $rO_jobClustering2 = RRNAAmplicons::clustering2(
				\%cfg,
				$fastq,
				$barcodes,
				$dir."/obs"
			);
			if(!$rO_jobClustering2->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "clustering", "clustering2", "CLUSTERING2"."_".$loop_counter, $dependency, "global", $rO_jobClustering2) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobClustering2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." clustering_method_2\n";
			$currStep++;

		}elsif($clusteringMethod == 3){
			my $rO_jobClustering3 = RRNAAmplicons::clustering3(
				\%cfg,
				$fastq,
				$barcodes,
				$dir."/obs"
			);
			if(!$rO_jobClustering3->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "clustering", "clustering3", "CLUSTERING3"."_".$loop_counter, $dependency, "global", $rO_jobClustering3) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobClustering3->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." clustering_method_3\n";
			$currStep++;
		}
		
		# RDP wrapper
		my $rO_jobRDPWrapper = MicrobialEcology::RDPWrapper(
			\%cfg,
			$dir."/obs/derep1_099_derep2_nonChimDeNovoRef_sorted_097_sorted.fasta",
			$dir."/rdp/rdp.tsv"
		);
		if(!$rO_jobRDPWrapper->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "RDP", "RDP", "RDP"."_".$loop_counter, $dependency, "global", $rO_jobRDPWrapper) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobRDPWrapper->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." RDP\n";
		$currStep++;
		
		# Add rdp taxonomy to SeqObs table
		my $rO_jobAddTaxToObs = MicrobialEcology::addTaxToObs(
			\%cfg,
			$dir."/obs/obs_filtered.tsv",
			$dir."/rdp/rdp.tsv",
			$dir."/otu_tables/otu_table.tsv",
			$dir."/otu_tables/otu_table_failed.tsv"
		);
		if(!$rO_jobAddTaxToObs->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "add_taxonomy", "add_taxonomy", "ADDTAXONOMY"."_".$loop_counter, $dependency, "global", $rO_jobAddTaxToObs) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobAddTaxToObs->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." add_taxonomy\n";
		$currStep++;
	
		# Filter OTU table (remove columns that contains all zeros). Will crash on upcoming steps if we don't do so here...
		my $rO_jobFilterOTUTable = MicrobialEcology::filterOTUTable(
			\%cfg,
			$dir."/otu_tables/otu_table.tsv",
			$dir."/otu_tables/otu_table_filtered.tsv" # filtered here means only that columns with zeros are being removed.
		);
		if(!$rO_jobFilterOTUTable->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "filter_obs", "filter_obs", "FILTEROBS"."_".$loop_counter, $dependency, "global", $rO_jobFilterOTUTable) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobFilterOTUTable->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." filter_obs\n";
		$currStep++;
	
			
		# Generate OTU table containing bacteria only or Fungi only depending...
		if($fungi_ITS){
			my $rO_jobSplitOTUTable = MicrobialEcology::splitOTUTable(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered.tsv",
				$dir."/otu_tables/otu_table_filtered_fungi.tsv",
				$dir."/otu_tables/otu_table_filtered_others.tsv",
				"fungi"
			);
			if(!$rO_jobSplitOTUTable->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "split_otu_table", "split_otu_table", "SPLITOTUTABLE"."_".$loop_counter, $dependency, "global", $rO_jobSplitOTUTable) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobSplitOTUTable->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." split_otu_table\n";
			$currStep++;

		}elsif($bactArch){ # bacteria/arch
			my $rO_jobSplitOTUTable = MicrobialEcology::splitOTUTable(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered.tsv",
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea.tsv",
				$dir."/otu_tables/otu_table_filtered_others.tsv",
				"bactArch"
			);
			if(!$rO_jobSplitOTUTable->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "split_otu_table", "split_otu_table", "SPLITOTUTABLE"."_".$loop_counter, $dependency, "global", $rO_jobSplitOTUTable) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobSplitOTUTable->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." split_otu_table\n";
			$currStep++;

		}elsif($everything){
			# Do nothing.

		}
	
		##################################
		# Convert OTU table to biom table
		#
		my ($tsv, $biom);
		if($fungi_ITS){
			$tsv = $dir."/otu_tables/otu_table_filtered_fungi.tsv";
			$biom = $dir."/otu_tables/otu_table_filtered_fungi.biom";
		}elsif($bactArch){
			$tsv = $dir."/otu_tables/otu_table_filtered_bacteriaArchaea.tsv";
			$biom = $dir."/otu_tables/otu_table_filtered_bacteriaArchaea.biom";
		}elsif($everything){

		}

		unless($everything){
			my $rO_jobConvertOTUToBiom = MicrobialEcology::convertOTUToBiom(
				\%cfg,
				$tsv,
				$biom
			);
			if(!$rO_jobConvertOTUToBiom->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "convert_otu_to_biom", "convert_otu_to_biom", "CONVERTOTUTOBIOM"."_".$loop_counter, $dependency, "global", $rO_jobConvertOTUToBiom) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobConvertOTUToBiom->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." convert_otu_to_biom\n";
			$currStep++;
		
			# Convert OTU table to biom table
			my $rO_jobConvertOTUToBiom2 = MicrobialEcology::convertOTUToBiom(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered_others.tsv",
				$dir."/otu_tables/otu_table_filtered_others.biom"
			);
				if(!$rO_jobConvertOTUToBiom2->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "convert_otu_to_biom", "convert_otu_to_biom2", "CONVERTOTUTOBIOM2"."_".$loop_counter, $dependency, "global", $rO_jobConvertOTUToBiom2) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobConvertOTUToBiom2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." convert_otu_to_biom2\n";
			$currStep++;
		}	
		
		# Convert OTU table to biom table raw otu table.
		my $rO_jobConvertOTUToBiom3 = MicrobialEcology::convertOTUToBiom(
			\%cfg,
			$dir."/otu_tables/otu_table_filtered.tsv",
			$dir."/otu_tables/otu_table_filtered.biom"
		);
		if(!$rO_jobConvertOTUToBiom3->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "convert_otu_to_biom", "convert_otu_to_biom3", "CONVERTOTUTOBIOM3"."_".$loop_counter, $dependency, "global", $rO_jobConvertOTUToBiom3) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobConvertOTUToBiom3->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." convert_otu_to_biom3\n";
		$currStep++;
		
		#####################################################################
		# Rarefy otu tables and convert resulting .biom files in .tsv files
		#
		if($fungi_ITS){
			my $rO_jobRarefy = MicrobialEcology::rarefy(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered_fungi.biom",
				$dir."/otu_tables/otu_table_filtered_fungi_rarefied.biom",
				$dir."/otu_tables/otu_table_filtered_fungi_rarefied.tsv"
			);
			if(!$rO_jobRarefy->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "rarefy", "rarefy", "RAREFY"."_".$loop_counter, $dependency, "global", $rO_jobRarefy) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobRarefy->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." rarefy\n";
			$currStep++;

		}elsif($bactArch){
			my $rO_jobRarefy = MicrobialEcology::rarefy(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea.biom",
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied.biom",
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied.tsv"
			);
			if(!$rO_jobRarefy->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "rarefy", "rarefy", "RAREFY"."_".$loop_counter, $dependency, "global", $rO_jobRarefy) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobRarefy->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." rarefy\n";
			$currStep++;
		}
	
		my $rO_jobRarefy = MicrobialEcology::rarefy(
			\%cfg,
			$dir."/otu_tables/otu_table_filtered.biom",
			$dir."/otu_tables/otu_table_filtered_rarefied.biom",
			$dir."/otu_tables/otu_table_filtered_rarefied.tsv"
		);
		if(!$rO_jobRarefy->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "rarefy", "rarefy", "RAREFY"."_".$loop_counter, $dependency, "global", $rO_jobRarefy) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobRarefy->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." rarefy\n";
		$currStep++;
		
		####################################################################################################
		# Generate a rarefied filtered otu table of 2x2. #TODO merge in a loop instead of duplicating code.
		#
		if($fungi_ITS){
			my $rO_jobFilter2X2 = MicrobialEcology::filterAndConvertToBiom(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered_fungi_rarefied.tsv",
				$dir."/otu_tables/otu_table_filtered_fungi_rarefied_2X2.tsv",
				2,
				2,
				$dir."/otu_tables/otu_table_filtered_fungi_rarefied_2X2.biom"
			);
			if(!$rO_jobFilter2X2->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "filter_obs_table", "filter_obs_table", "FILTEROBSTABLE"."_".$loop_counter, $dependency, "global", $rO_jobFilter2X2) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobFilter2X2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." filter_obs_table\n";
			$currStep++;

		}elsif($bactArch){
			my $rO_jobFilter2X2 = MicrobialEcology::filterAndConvertToBiom(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied.tsv",
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied_2X2.tsv",
				2,
				2,
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied_2X2.biom"
			);
			if(!$rO_jobFilter2X2->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "filter_obs_table", "filter_obs_table2X2", "FILTEROBSTABLE2X2"."_".$loop_counter, $dependency, "global", $rO_jobFilter2X2) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobFilter2X2->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." filter_obs_table2X2\n";
			$currStep++;
		}
		
		my $rO_jobFilter2X2All = MicrobialEcology::filterAndConvertToBiom(
			\%cfg,
			$dir."/otu_tables/otu_table_filtered_rarefied.tsv",
			$dir."/otu_tables/otu_table_filtered_rarefied_2X2.tsv",
			2,
			2,
			$dir."/otu_tables/otu_table_filtered_rarefied_2X2.biom"
		);
		if(!$rO_jobFilter2X2All->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "filter_obs_table", "filter_obs_table2X2_all", "FILTEROBSTABLE2X2ALL"."_".$loop_counter, $dependency, "global", $rO_jobFilter2X2All) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobFilter2X2All->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." filter_obs_table2X2\n";
		$currStep++;

		# Now the 1X1 tables	
		if($fungi_ITS){
			my $rO_jobFilter1X1 = MicrobialEcology::filterAndConvertToBiom(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered_fungi_rarefied.tsv",
				$dir."/otu_tables/otu_table_filtered_fungi_rarefied_1X1.tsv",
				1,
				1,
				$dir."/otu_tables/otu_table_filtered_fungi_rarefied_1X1.biom"
			);
			if(!$rO_jobFilter1X1->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "filter_obs_table", "filter_obs_table1X1", "FILTEROBSTABLE1X1"."_".$loop_counter, $dependency, "global", $rO_jobFilter1X1) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobFilter1X1->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." filter_obs_table1X1\n";
			$currStep++;
				
		}elsif($bactArch){
			my $rO_jobFilter1X1 = MicrobialEcology::filterAndConvertToBiom(
				\%cfg,
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied.tsv",
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied_1X1.tsv",
				1,
				1,
				$dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied_1X1.biom"
			);
			if(!$rO_jobFilter1X1->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "filter_obs_table", "filter_obs_table1X1", "FILTEROBSTABLE1X1"."_".$loop_counter, $dependency, "global", $rO_jobFilter1X1) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobFilter1X1->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." filter_obs_table1X1\n";
			$currStep++;
		}
		
		my $rO_jobFilter1X1All = MicrobialEcology::filterAndConvertToBiom(
			\%cfg,
			$dir."/otu_tables/otu_table_filtered_rarefied.tsv",
			$dir."/otu_tables/otu_table_filtered_rarefied_1X1.tsv",
			1,
			1,
			$dir."/otu_tables/otu_table_filtered_rarefied_1X1.biom"
		);
		if(!$rO_jobFilter1X1All->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "filter_obs_table", "filter_obs_table1X1_all", "FILTEROBSTABLE1X1ALL"."_".$loop_counter, $dependency, "global", $rO_jobFilter1X1All) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobFilter1X1All->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." filter_obs_table1X1\n";
		$currStep++;

		#################################################################	
		# Define what OTU table to use for the remaining of the pipeline.
		# But also plot raw OTU table.
		#
		my $otu_table;
		my $otu_table_prefix;
		if($fungi_ITS){
			$otu_table = $dir."/otu_tables/otu_table_filtered_fungi_rarefied.biom";
			$otu_table_prefix = "otu_table_filtered_fungi_rarefied";
		}elsif($bactArch){
			$otu_table = $dir."/otu_tables/otu_table_filtered_bacteriaArchaea_rarefied.biom";	
			$otu_table_prefix = "otu_table_filtered_bacteriaArchaea_rarefied";
		}elsif($everything){
			$otu_table = $dir."/otu_tables/otu_table_filtered_rarefied.biom";	
			$otu_table_prefix = "otu_table_filtered_rarefied";
		}
		my $otu_table_raw = $dir."/otu_tables/otu_table_filtered.biom";
		my $otu_table_prefix_raw = "otu_table_filtered";
	
		# Summarize taxonomy with absolute abundance
    my $taxDepth = LoadConfig::getParam(\%cfg, 'summarize_taxonomy', 'taxonomyDepth', 1, 'int'); 
    my $plotTaxaStringAbs = "";
		for(my $i=1; $i<$taxDepth; $i++){
			my $rO_jobSummarizeTaxonomyAbsolute = MicrobialEcology::summarizeTaxonomyAbsolute(
				\%cfg,
				$otu_table,
				$i,
				$dir."/tax_summary/absolute/"
			);
			if(!$rO_jobSummarizeTaxonomyAbsolute->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "summarize_taxonomy", "summarize_taxonomy_absolute", "SUMMARIZETAXONOMYABSOLUTE_L$i"."_".$loop_counter, $dependency, "global", $rO_jobSummarizeTaxonomyAbsolute) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobSummarizeTaxonomyAbsolute->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." summarize_taxonomy_absolute_L".$i."\n";
      
      # construct string for next step
      $plotTaxaStringAbs .= $dir."/tax_summary/absolute/".$otu_table_prefix."_L".$i.".txt,"; 
			$currStep++;
		}
    chop($plotTaxaStringAbs);
  	
		# Summarize taxonomy with relative abundance
    my $plotTaxaStringRel = "";
		for(my $i=1; $i<$taxDepth; $i++){
			my $rO_jobSummarizeTaxonomyRelative = MicrobialEcology::summarizeTaxonomyRelative(
				\%cfg,
				$otu_table,
				$i,
				$dir."/tax_summary/relative/"
			);
			if(!$rO_jobSummarizeTaxonomyRelative->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "summarize_taxonomy", "summarize_taxonomy_relative", "SUMMARIZETAXONOMYRELATIVE_L$i"."_".$loop_counter, $dependency, "global", $rO_jobSummarizeTaxonomyRelative) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobSummarizeTaxonomyRelative->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." summarize_taxonomy_relative_L".$i."\n";
      
      # construct string for next step
      $plotTaxaStringRel .= $dir."/tax_summary/relative/".$otu_table_prefix."_L".$i.".txt,"; 
			$currStep++;
		}
    chop($plotTaxaStringRel);
		
		# Summarize taxonomy with absolute abundance - raw
    my $plotTaxaStringAbsRaw = "";
		for(my $i=1; $i<$taxDepth; $i++){
			my $rO_jobSummarizeTaxonomyAbsolute = MicrobialEcology::summarizeTaxonomyAbsolute(
				\%cfg,
				$otu_table_raw,
				$i,
				$dir."/tax_summary/absolute/"
			);
			if(!$rO_jobSummarizeTaxonomyAbsolute->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "summarize_taxonomy", "summarize_taxonomy_absolute_raw", "SUMMARIZETAXONOMYABSOLUTERAW"."_".$loop_counter, $dependency, "global", $rO_jobSummarizeTaxonomyAbsolute) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobSummarizeTaxonomyAbsolute->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." summarize_taxonomy_absolute_raw_L".$i."\n";
      
      # construct string for next step
      $plotTaxaStringAbsRaw .= $dir."/tax_summary/absolute/".$otu_table_prefix_raw."_L".$i.".txt,"; 
			
      $currStep++;
		}
    chop($plotTaxaStringAbsRaw);
		
		# Summarize taxonomy with relative abundance - raw
    my $plotTaxaStringRelRaw = "";
		for(my $i=1; $i<$taxDepth; $i++){
			my $rO_jobSummarizeTaxonomyRelative = MicrobialEcology::summarizeTaxonomyRelative(
				\%cfg,
				$otu_table_raw,
				$i,
				$dir."/tax_summary/relative/"
			);
			if(!$rO_jobSummarizeTaxonomyRelative->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "summarize_taxonomy", "summarize_taxonomy_relative_raw", "SUMMARIZETAXONOMYRELATIVERAW"."_".$loop_counter, $dependency, "global", $rO_jobSummarizeTaxonomyRelative) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobSummarizeTaxonomyRelative->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." summarize_taxonomy_relative_raw_L".$i."\n";

      # construct string for next step
      $plotTaxaStringRelRaw .= $dir."/tax_summary/relative/".$otu_table_prefix_raw."_L".$i.".txt,"; 

			$currStep++;
		}
    chop($plotTaxaStringRelRaw);
	
		#####################################################	
		# Qiime make taxa plot absolute abundance
		# Do it for both bacteria/archaea/fungi and raw table
		#
		my $rO_jobPlotTaxa = MicrobialEcology::plotTaxa(
			\%cfg,
			#$dir."/tax_summary/absolute/".$otu_table_prefix."_L1.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix."_L2.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix."_L3.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix."_L4.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix."_L5.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix."_L6.txt",
			$plotTaxaStringAbs,
			$dir."/tax_summary/plots"	
		);
		if(!$rO_jobPlotTaxa->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "plot_taxonomy", "plot_taxonomy_absolute", "PLOTTAXONOMYABSOLUTE"."_".$loop_counter, $dependency, "global", $rO_jobPlotTaxa) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobPlotTaxa->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." plot_taxonomy_absolute\n";
		$currStep++;
		
		my $rO_jobPlotTaxaRaw = MicrobialEcology::plotTaxa(
			\%cfg,
			#$dir."/tax_summary/absolute/".$otu_table_prefix_raw."_L1.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix_raw."_L2.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix_raw."_L3.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix_raw."_L4.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix_raw."_L5.txt,".$dir."/tax_summary/absolute/".$otu_table_prefix_raw."_L6.txt",
      $plotTaxaStringAbsRaw,
			$dir."/tax_summary/plots"	
		);
		if(!$rO_jobPlotTaxaRaw->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "plot_taxonomy", "plot_taxonomy_absolute_raw", "PLOTTAXONOMYABSOLUTERAW"."_".$loop_counter, $dependency, "global", $rO_jobPlotTaxaRaw) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobPlotTaxaRaw->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." plot_taxonomy_absolute_raw\n";
		$currStep++;

		# Calculate abundance thresholds
		#my $rO_jobCalculateAbundanceThresholds = MicrobialEcology::calculateAbundanceThresholds(
		#	\%cfg,
		#	$dir."/otu_tables/otu_table_filtered.tsv",
		#	$dir."/otu_tables/"
		#);
		#if(!$rO_jobCalculateAbundanceThresholds->isUp2Date()) {
		#	SubmitToCluster::printSubmitCmd(\%cfg, "abundance_thresholds", "abundance_thresholds_absolute", "ABUNDANCETHRESHOLDS"."_".$loop_counter, $dependency, "global", $rO_jobCalculateAbundanceThresholds) if($start_at <= $currStep && $end_at >= $currStep); 
		#	$dependency = $rO_jobCalculateAbundanceThresholds->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		#}
		#print STDERR "[DEBUG]\t\tStep #".$currStep."\n";
		#$currStep++;
	
		# Log rank abundance
		# Remove $dir."log_rank_abundance dir if exists.
		#$cmd = 'module load '.LoadConfig::getParam(\%cfg, 'memtime', 'moduleVersion.memtime').' ; memtime ';	 
		#$cmd .= $plotRankAbundance." ";
		#$cmd .= "-s '*' ";
		#$cmd .= "-i ".$dir."/otu_tables/otu_table_filtered.biom ";
		#$cmd .= "-o ".$dir."/log_rank_abundance ";
		#if(-d $dir."/log_rank_abundance"){
		#	system("rm -rf ".$dir."/log_rank_abundance");
		#}
		#$dependency = executeCmd($cmd, \%cfg, "log_rank_abundance", "Log_rank_abundance_".$loop_counter, $dependency, $currStep) if($start_at <= $currStep && $end_at >= $currStep);
		#print STDERR "[DEBUG]\t\tStep #".$currStep."\n";
		#$currStep++;
		
		# generate a phylum barplot for a quick view of how the data looks like.
		my $rO_jobPhylumBarplot = MicrobialEcology::phylumBarplot(
			\%cfg,
			$dir."/tax_summary/absolute/otu_table_filtered_L2.txt",
			$dir."/tax_summary/",
			"taxonomy_phylum_L2"
		);
		if(!$rO_jobPhylumBarplot->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "phylum_barplot", "phylum_barplot_all", "PHYLUMBARPLOTALL"."_".$loop_counter, $dependency, "global", $rO_jobPhylumBarplot) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobPhylumBarplot->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." phylum_barplot_all\n";
		$currStep++;
		
		my $rO_jobPhylumBarplotRel = MicrobialEcology::phylumBarplot(
			\%cfg,
			$dir."/tax_summary/relative/otu_table_filtered_L2.txt",
			$dir."/tax_summary/",
			"taxonomy_phylum_L2_relative"
		);
		if(!$rO_jobPhylumBarplotRel->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "phylum_barplot", "phylum_barplot_all_rel", "PHYLUMBARPLOTREL"."_".$loop_counter, $dependency, "global", $rO_jobPhylumBarplotRel) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobPhylumBarplotRel->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\t\tStep #".$currStep." phylum_barplot_all_rel\n";
		$currStep++;
	
		# generate a phylum barplot for a quick view of how the data looks like (bacteria and archaea only).
		if($fungi_ITS){	
			my $rO_jobPhylumBarplotFungi = MicrobialEcology::phylumBarplot(
				\%cfg,
				$dir."/tax_summary/absolute/otu_table_filtered_fungi_rarefied_L2.txt",
				$dir."/tax_summary/",
				"taxonomy_phylum_filtered_fungi_rarefied_L2"
			);
			if(!$rO_jobPhylumBarplotFungi->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "phylum_barplot", "phylum_barplot_fungi_rarefied", "PHYLUMBARPLOTFUNGI"."_".$loop_counter, $dependency, "global", $rO_jobPhylumBarplotFungi) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobPhylumBarplotFungi->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." phylum_barplot_fungi_rarefied\n";
			$currStep++;

		}elsif($bactArch){
			my $rO_jobPhylumBarplotBacteria = MicrobialEcology::phylumBarplot(
				\%cfg,
				$dir."/tax_summary/absolute/otu_table_filtered_bacteriaArchaea_rarefied_L2.txt",
				$dir."/tax_summary/",
				"taxonomy_phylum_bacteria_archaea_L2"
			);
			if(!$rO_jobPhylumBarplotBacteria->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "phylum_barplot", "phylum_barplot_bacteria_rarefied", "PHYLUMBARPLOTBACTERIA"."_".$loop_counter, $dependency, "global", $rO_jobPhylumBarplotBacteria) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobPhylumBarplotBacteria->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." phylum_barplot_bacteria_rarefied\n";
			$currStep++;

		}elsif($everything){
			my $rO_jobPhylumBarplotAll = MicrobialEcology::phylumBarplot(
				\%cfg,
				$dir."/tax_summary/absolute/otu_table_filtered_rarefied_L2.txt",
				$dir."/tax_summary/",
				"taxonomy_phylum_all_L2"
			);
			if(!$rO_jobPhylumBarplotAll->isUp2Date()) {
				SubmitToCluster::printSubmitCmd(\%cfg, "phylum_barplot", "phylum_barplot_all_rarefied", "PHYLUMBarplotAll"."_".$loop_counter, $dependency, "global", $rO_jobPhylumBarplotAll) if($start_at <= $currStep && $end_at >= $currStep); 
				$dependency = $rO_jobPhylumBarplotAll->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
			}
			print STDERR "[DEBUG]\t\tStep #".$currStep." phylum_barplot_all_rarefied\n";
			$currStep++;
		}

    ###############################
    # Blastn - raw OTUs. If PacBio only
    #
    my $dataType = LoadConfig::getParam(\%cfg, 'default', 'sequencer');
    if($dataType =~ m/pacbio/i){
	    my $rO_jobBlast = MicrobialEcology::blast(
  		  \%cfg,
  			$dir."/obs/obs_filtered.fasta",
  			$dir."/blast/obs_filtered.blastout"
  		);
  		if(!$rO_jobBlast->isUp2Date()) {
  			SubmitToCluster::printSubmitCmd(\%cfg, "blast", "blast_raw_otus", "BLASTRAWOTUS"."_".$loop_counter, $dependency, "global", $rO_jobBlast) if($start_at <= $currStep && $end_at >= $currStep); 
  			$dependency = $rO_jobBlast->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
  		}
  		print STDERR "[DEBUG]\t\tStep #".$currStep." blast raw OTUs\n";
  		$currStep++;
	  }	

		return $dependency;
	}
		
	# sub makedirs: Make directories for each type of analyses. 
	# @input: a directory
	# @returns: nothin'!
	sub make_dirs{
		my $prefix = $_[0];	
		
		print STDOUT "mkdir -p ".$prefix."\n";	
		system("mkdir -p ".$prefix) if $noMsub;	

		if($prefix =~ m/\/assembled/){
			#print STDOUT "mkdir -p ".$prefix."/fastqs\n";
			#system("mkdir -p ".$prefix."/fastqs") if $noMsub;
			
		}
		my $qscore_dir = $prefix."/qscores";
		print STDOUT "mkdir -p $qscore_dir\n";
		system("mkdir -p ".$qscore_dir) if $noMsub;

		my $obs_dir = $prefix."/obs\n";
		print STDOUT "mkdir -p ".$obs_dir;
		system("mkdir -p ".$obs_dir) if $noMsub;

		my $rdp_dir = $prefix."/rdp\n";
		print STDOUT "mkdir -p ".$rdp_dir;
		system("mkdir -p ".$rdp_dir) if $noMsub;
	
		my $otu_table_dir = $prefix."/otu_tables\n";
		print STDOUT "mkdir -p ".$otu_table_dir;
		system("mkdir -p ".$otu_table_dir) if $noMsub;
		
    my $blast_dir = $prefix."/blast\n";
		print STDOUT "mkdir -p ".$blast_dir;
		system("mkdir -p ".$blast_dir) if $noMsub;

		my $tax_dir = $prefix."/tax_summary";
		print STDOUT "mkdir -p ".$tax_dir."\n";
		print STDOUT "mkdir -p ".$tax_dir."/absolute\n";
		print STDOUT "mkdir -p ".$tax_dir."/relative\n";
		print STDOUT "mkdir -p ".$tax_dir."/plots\n";
		system("mkdir -p ".$tax_dir) if $noMsub;
		system("mkdir -p ".$tax_dir."/absolute") if $noMsub;
		system("mkdir -p ".$tax_dir."/relative") if $noMsub;
		system("mkdir -p ".$tax_dir."/plots") if $noMsub;
	}

	# Merge pdf
	# @input: files to merges
	# $output: pdf file
			#$curr_dir."/countReport.pdf",
			#$qscore_pdf,
			#$curr_dir."/REPORT_".$barcodes_filename.".pdf",
			#"ex"
	sub mergePdf{
		my $dependency = shift;
		my $prefix = shift(@_);
		my $file_1 = shift(@_);
		my $outfile = pop(@_);
		
		my $cmd = "";
		$cmd .= "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=".$outfile;
		$cmd .= " ".$file_1;

		foreach my $file (@_){
			$cmd .= " ".$file;
		}
	
		my $rO_jobmergePdf = RRNAAmplicons::mergePdf(
			\%cfg,
			$cmd
		);
		if(!$rO_jobmergePdf->isUp2Date()) {
			SubmitToCluster::printSubmitCmd(\%cfg, "merge_pdf", "merge_pdf", "MERGEPDF"."_".$loop_counter, $dependency, "global", $rO_jobmergePdf) if($start_at <= $currStep && $end_at >= $currStep); 
			$dependency = $rO_jobmergePdf->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
		}
		print STDERR "[DEBUG]\tStep #".$currStep." merge_pdf\n";
		$currStep++;
		return $dependency;
	}

	# =========================== #
	# ==== SUBROUTINES END ====== #	
	# =========================== #
}

# Set script name (without suffix) as pipeline name
my $pipelineName = fileparse($0, qr/\.[^.]*/) . "-$Version::version";
my $stepNames = join(",", map($steps[$_]->{'name'}, @stepRange));

# Log anynymous statistics on remote MUGQIC web server
Tools::mugqicLog($pipelineName, $stepNames, $currNumberOfSamples);

exit;

# remove temp files
sub END{

  if(defined $currStep){
  	local $?;
  	my $rO_jobCleanup = RRNAAmplicons::cleanup(
  		\%cfg,
  		$tmpdir
  	);
  	if(!$rO_jobCleanup->isUp2Date()) {
  		SubmitToCluster::printSubmitCmd(\%cfg, "cleanup", "cleanup", "CLEANUP", $dependency, "global", $rO_jobCleanup) if($start_at <= $currStep && $end_at >= $currStep); 
  		$dependency = $rO_jobCleanup->getCommandJobId(0) if($start_at <= $currStep && $end_at >= $currStep);
  	}
  	print STDERR "[DEBUG]\tStep #".$currStep." cleanup\n";
  	$currStep++;
  }
}

