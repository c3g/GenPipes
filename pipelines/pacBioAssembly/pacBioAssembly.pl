#!/usr/bin/env perl

=head1 NAME

I<pacBioAssembly>

=head1 SYNOPSIS

pacBioAssembly.pl

=head1 DESCRIPTION

B<rnaSeq> Is the main PacBio HGAP assembly.

=head1 AUTHOR

B<Julien Tremblay> - I<julien.tremblay@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing

=cut

# Strict Pragmas
#---------------------
use strict;
use warnings;
#---------------------

BEGIN{
    #Makesure we can find the GetConfig::LoadModules module relative to this script install
    use File::Basename;
    use Cwd 'abs_path';
    my ( undef, $mod_path, undef ) = fileparse( abs_path(__FILE__) );
    unshift @INC, $mod_path."lib";
	#$ENV{PATH} = $ENV{PATH}.":".$mod_path;
	#print STDERR $ENV{PATH}."\n";
	#print STDERR join(":", @INC); print STDERR "\n";
}


# Dependencies
#--------------------
use Getopt::Std;
use Cwd qw/ abs_path /;

use LoadConfig;
use SmrtAnalysis;
use PrinSeq;
use PacBioTools;
use Celera;
use SampleSheet;
use SequenceDictionaryParser;
use SubmitToCluster;
use BLAST;
use Mummer;
use GqSeqUtils;
#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'filtering', 	'stepLoop' => 'global',    'output' => 'filtering'});
push(@steps, {'name' => 'getStats', 	'stepLoop' => 'sample',    'output' => 'getStats'});
push(@steps, {'name' => 'preassembly', 	'stepLoop' => 'sample',    'output' => 'preassembly'});
push(@steps, {'name' => 'assembly', 	'stepLoop' => 'assembly',  'output' => 'assembly'});
push(@steps, {'name' => 'polishing',	'stepLoop' => 'assembly',  'output' => 'realign'});
push(@steps, {'name' => 'blast',        'stepLoop' => 'assembly',  'output' => 'blastContigs'});
push(@steps, {'name' => 'mummer', 		'stepLoop' => 'assembly',  'output' => 'mummer'});
push(@steps, {'name' => 'report', 		'stepLoop' => 'assembly',  'output' => 'report'});

my %globalDep;
for my $stepName (@steps) { 
	$globalDep{$stepName -> {'name'} } ={};
}

# Global scope variables
my $designFilePath;
my $configFile;
my $workDir;
my $smrtCells;
my $samplesFile;
my $noMsub;
my $estimatedGenomeSize;
my $numberOfBases;
my $libType;
my $verbose = 0;
my $TMPDIR;
my $hgapAlgorithm;

&main();

sub printUsage {
  print "\nUsage: perl ".$0." \n";
  print "\t-c  config file\n";
  print "\t-s  start step, inclusive\n";
  print "\t-e  end step, inclusive\n";
  print "\t-n  nanuq sample sheet\n";
  print "\t-w  work directory\n";
  print "\t-m  flag to run job interactively\n";
  print "\t-h  1: hgap1 legacy(smrtpipe). 2: hgap1. 3:hgap2\n";
  print "\t-v  verbose\n";
  print "\n";
  print "Steps:\n";
  for(my $idx = 0; $idx < @steps; $idx++) {
    print "".($idx + 1) . '- ' . $steps[$idx]->{'name'} . "\n";
  }
  print "\n";
}

sub main {
	my %opts;
	getopts('c:s:e:n:w:m:v:h:', \%opts);
	
	if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'e'}) || !defined($opts{'n'}) || !defined($opts{'w'} ) ) {
		printUsage();
		exit(1);
	}
	
	# Load config file
	$noMsub = 1 if(defined($opts{'m'}));
	my %jobIdVarPrefix;
	my %cfg 				= LoadConfig->readConfigFile($opts{'c'});
	$TMPDIR				    = LoadConfig::getParam(\%cfg, 'default', 'tmpDir');	

	# Get design groups
	$workDir 				= abs_path($opts{'w'});
	$configFile 			= abs_path($opts{'c'});
	$samplesFile			= abs_path($opts{'n'}); 
	$hgapAlgorithm          = $opts{'h'};
	$hgapAlgorithm          = 0 if(!defined($hgapAlgorithm));
	# TODO add a step to validate values found in config file.

	# Define outdir
	my $outdir = LoadConfig::getParam(\%cfg, 'default', 'outdir');

	# Mkdir for fofns
	system("mkdir -p $outdir/fofns") ;#if$noMsub;	
	print STDOUT "mkdir -p $outdir/fofns\n";

	# Loop through samples. Could only use one hash...
	my %hohSamples;
	my %hEstimatedGenomeSize;
	my %hSmrtCells;
	my %hLibType;
	my %hNumberOfBases;

	open(SAMPLES, '<'.$samplesFile) or die "Can't find file $samplesFile\n";
	while(<SAMPLES>){
		chomp;
		my @row = split(/\t/, $_);
		$hohSamples{$row[0]}{$row[1]}  = $row[2];	
		$hLibType{$row[0]} = $row[3];
		if(exists $hNumberOfBases{$row[0]}){
			$hNumberOfBases{$row[0]} = $hNumberOfBases{$row[0]} + $row[4];
		}else{
			$hNumberOfBases{$row[0]} = $row[4];	
		}
		$hEstimatedGenomeSize{$row[0]} = $row[5];
		$hSmrtCells{$row[0]}++;
	}
	close(SAMPLES);	
	
	# fofns will be written here. One fofn per sample. One smrt cell can the same sample more than one time.
	my %fofns;
	while(my($sampleName,$value1)=each(%hohSamples)){
		#print STDERR "$sampleName\n";		
		
		# Write fofns here in Main
		open(FOFN, ">$outdir/fofns/$sampleName.fofn") or die "Can't open fofn file for writing ... $!\n";
		
		while(my($well,$pathToFofn)=each(%$value1)){
			#print STDERR $sampleName."\t".$well."\t".$pathToFofn."\n";
			print FOFN "$pathToFofn\n";
  		}
		close(FOFN);

		$fofns{$sampleName} = ">$outdir/fofns/$sampleName.fofn";
	}

	SubmitToCluster::initPipeline;

	for my $currSampleName ( keys %fofns ) {
		my $currFofn = $fofns{$currSampleName};
    		
		# Define estimated genome Size in global variable.
		$estimatedGenomeSize = $hEstimatedGenomeSize{$currSampleName};
		$smrtCells           = $hSmrtCells{$currSampleName};
		$libType             = $hLibType{$currSampleName};		
		$numberOfBases       = $hNumberOfBases{$currSampleName};		
	
		# Find coverage range using estmated genome size and number of bases.
		my $estimatedCoverage = int($numberOfBases / $estimatedGenomeSize);
		
		print STDERR "[DEBUG] EstimatedGenomeSize: ".$estimatedGenomeSize."\n";
		print STDERR "[DEBUG] smrtCells: ".$smrtCells."\n";
		print STDERR "[DEBUG] LibType: ".$libType."\n";
		print STDERR "[DEBUG] numberOfbases: ".$numberOfBases."\n";
		print STDERR "[DEBUG] Estimated coverage: ".$estimatedCoverage."\n";	

		my $currentStep;
		my $rootDependency = undef;
		my $dependency = undef;
		# Call core pipeline methods here.
		for($currentStep = $opts{'s'}-1; $currentStep <= ($opts{'e'}-1); $currentStep++) {
			#$rootDependency = filtering("filtering", $currSampleName, $currFofn, \%cfg, $dependency);

  			my $fname = $steps[$currentStep]->{'name'};
   			my $subref = \&$fname;
     
			if ($steps[$currentStep]->{'stepLoop'} eq 'global') {
				# Tests for the first step in the list. Used for dependencies.
				$rootDependency = &$subref($steps[$currentStep]->{'name'}, $currSampleName, $currFofn, \%cfg, $dependency); 
				$globalDep{$fname}->{$currSampleName} = $rootDependency;
			}   
		}   
		
		# do assembly for the various coverage values entered in the ini file. And also do assembly with
		# the "logical" total sequenced bp / genome estimated size.	
		my $coverageRange = LoadConfig::getParam(\%cfg, 'stats', 'estimatedCoverage');
		my @coverageRange = split(/_/, $coverageRange);
		push(@coverageRange, $estimatedCoverage);
		
		my $merSizes = LoadConfig::getParam(\%cfg, 'default', 'merSizes');
		my @merSizes = split(/:/, $merSizes);

		# Sample
		foreach my $coverage(@coverageRange){
	
			for($currentStep = $opts{'s'}-1; $currentStep <= ($opts{'e'}-1); $currentStep++) {
	
  				my $fname = $steps[$currentStep]->{'name'};
   				my $subref = \&$fname;
   	  
				if ($steps[$currentStep]->{'stepLoop'} eq 'sample') {
					# Tests for the first step in the list. Used for dependencies.
					$dependency = &$subref($steps[$currentStep]->{'name'}, $currSampleName, $coverage."X", $currFofn, \%cfg, $rootDependency); 
					$globalDep{$fname}->{$currSampleName} = $dependency;

				}elsif($steps[$currentStep]->{'stepLoop'} eq 'assembly') {
					foreach my $merSize(@merSizes){
						print STDERR "[DEBUG] MERSIZE: ".$merSize."\n";
						# Tests for the first step in the list. Used for dependencies.
						$dependency = &$subref($steps[$currentStep]->{'name'}, $currSampleName, $coverage."X", $merSize, $currFofn, \%cfg, $rootDependency); 
						$globalDep{$fname}->{$currSampleName} = $dependency;
					}
				} 
			}  	
		}
	}
	exit;		
}

############
# Filtering. The RS_Filter_Only protocol pretty much does only one thing: generate fastq and subreads. 
# It also generates informative run metrics such as loading efficiency, readlengths, and base quality.
# @input   path.fofn: file containing one line pointing to a .h5 file.
#          filter_settings.xml
# @output  input.xml, outfiles(fastqs) in samplename/data/, smrtpipe.log
sub filtering{
	
	my $stepName	= shift;
	my $sampleName 	= shift;
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;
	
	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print "mkdir -p $outdir/$sampleName\n";
	print "mkdir -p $outdir/$sampleName/filtering\n";
	system("mkdir -p $outdir/$sampleName") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/filtering") ;#if$noMsub;
	
	# create PacBio input.xml
	my $rO_jobFofn = SmrtAnalysis::fofns(
		$rH_cfg, 
		"$outdir/fofns/$sampleName.fofn", 
		"$outdir/$sampleName/filtering/input.xml",
		"$outdir/$sampleName/filtering/input.fofn"
	);
	if(!$rO_jobFofn->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "fofn", $stepName , 'FILTERING_FOFNS', $dependency, $sampleName, $rO_jobFofn); 
	}

	# smart pipe to create fastq in --output <dir> arg.
	my $rO_jobFiltering = SmrtAnalysis::run(
		$rH_cfg,
		LoadConfig::getParam($rH_cfg, 'default', 'filteringSettings'),
		"$outdir/$sampleName/filtering/input.xml",
		"$outdir/$sampleName/filtering",
		"$outdir/$sampleName/filtering/smrtpipe.log"
	);
	if(!$rO_jobFiltering->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "smrtanalysis", $stepName , 'FILTERING_SMRTPIPE', $rO_jobFofn->getCommandJobId(0), $sampleName, $rO_jobFiltering); 
	}

	# Fastq to Fasta
	my $rO_jobPrinSeq = PrinSeq::fastqToFasta(
		$rH_cfg,
		"$outdir/$sampleName/filtering/data/filtered_subreads.fastq",
		"$outdir/$sampleName/filtering/data/filtered_subreads"
	);
	if(!$rO_jobPrinSeq->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "fastqToFasta", $stepName , 'FILTERING_FASTQTOFASTA', $rO_jobFiltering->getCommandJobId(0), $sampleName, $rO_jobPrinSeq); 
	}
 
	return $rO_jobPrinSeq->getCommandJobId(0);	
}

#################
# getStats
#
sub getStats{
	
	my $stepName	= shift;
	my $sampleName 	= shift;
	my $suffix	 	= shift; #will be ?X coverage value.
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;
	
	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/preassembly\n";
	system("mkdir -p $outdir/$sampleName/$suffix") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/$suffix/preassembly") ;#if$noMsub;

	my $coverage = $suffix; $coverage =~ s/(\d+)X/$1/;	
	
	# Choose a subread length threshold such that subreads above the threshold provide about 20x coverage of the genome.
	my $rO_jobGetCutoff = PacBioTools::getCutoff(
		$rH_cfg,
		"$outdir/$sampleName/filtering/data/filtered_subreads.fasta",
		$coverage,
		#LoadConfig::getParam($rH_cfg, 'default', 'genomeSize'), #TODO Put genomeSize value in sample sheet.	
		$estimatedGenomeSize,
		LoadConfig::getParam($rH_cfg, 'default', 'preassemblySettings'),
		"$outdir/$sampleName/$suffix/preassembly.xml",
		"$outdir/$sampleName/$suffix/preassemblyMinReadSize.txt"
	);
	if(!$rO_jobGetCutoff->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "pacbio_tools", $stepName , 'GETCUTOFF', $dependency, $sampleName."_".$suffix, $rO_jobGetCutoff); 
	}
 
	return $rO_jobGetCutoff->getCommandJobId(0);
}

#################
# Preassembly 
# Will run P_filter protocol followed by P_PreAssembler which includes alignment with BLASR
#
sub preassembly{
	my $stepName	= shift;
	my $sampleName 	= shift;
	my $suffix	 	= shift; #will be ?X coverage value.
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;

	my $cmd;

	# create PacBio input.xml 	
	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/preassembly\n";
	system("mkdir -p $outdir/$sampleName/$suffix") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/$suffix/preassembly") ;#if$noMsub;

	if($hgapAlgorithm == 0){	
		# create PacBio input.xml
		my $rO_jobFofn = SmrtAnalysis::fofns(
			$rH_cfg,
			"$outdir/fofns/$sampleName.fofn",
			"$outdir/$sampleName/$suffix/preassembly/input.xml",
			"/dev/null"
		);
		if(!$rO_jobFofn->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "fofn", $stepName , 'FOFNS', $dependency, $sampleName."_".$suffix, $rO_jobFofn); 
		}
	
		# smart pipe to create fastq in --output <dir> arg.
		my $rO_jobPreassembly = SmrtAnalysis::run(
			$rH_cfg,
			#"$outdir/$sampleName/$suffix/preassembly.xml",
			LoadConfig::getParam($rH_cfg, 'default', 'preassemblySettings'),
			"$outdir/$sampleName/$suffix/preassembly/input.xml",
			"$outdir/$sampleName/$suffix/preassembly",
			"$outdir/$sampleName/$suffix/preassembly/smrtpipe.log"
		);
		if(!$rO_jobPreassembly->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "smrtanalysis", $stepName , 'PREASSEMBLY', $rO_jobFofn->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobPreassembly); 
		}
	
		# Fastq to Fasta
		my $rO_jobPrinSeq = PrinSeq::fastqToFasta(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/preassembly/data/corrected.fastq",
			"$outdir/$sampleName/$suffix/preassembly/data/corrected"
		);
		if(!$rO_jobPrinSeq->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "fastqToFasta", $stepName, 'FASTQTOFASTA', $rO_jobPreassembly->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobPrinSeq); 
		}
	
		return $rO_jobPrinSeq->getCommandJobId(0);	
	
	}elsif($hgapAlgorithm == 1){

		# Separate Long reads from other reads.
		my $rO_jobSplitReads = PacBioTools::splitReads(
			$rH_cfg,
			"$outdir/$sampleName/filtering/data/filtered_subreads.fasta",
			"$outdir/$sampleName/$suffix/preassemblyMinReadSize.txt",	
			"$outdir/$sampleName/filtering/data/filtered_shortreads.fa",	
			"$outdir/$sampleName/filtering/data/filtered_longreads.fa"
		);
		if(!$rO_jobSplitReads->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "splitReads", $stepName, 'SPLITREADS', $dependency, $sampleName."_".$suffix, $rO_jobSplitReads); 
		}	
	
		# blasr
		my $rO_jobBlasr = SmrtAnalysis::blasr(
			$rH_cfg,
			"$outdir/$sampleName/filtering/data/filtered_subreads.fa",
			"$outdir/$sampleName/filtering/data/filtered_longreads.fa",
			"$outdir/$sampleName/$suffix/preassembly/seeds.sam",
			"$outdir/$sampleName/$suffix/preassembly/seeds.sam.fofn",
			"sam"
		);
		if(!$rO_jobBlasr->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "blasr", $stepName, 'BLASR', $rO_jobSplitReads->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobBlasr); 
		}
		
		# samtoh5
		my $rO_jobSamtoh5 = SmrtAnalysis::samtoh5(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/preassembly/seeds.sam",
			"$outdir/$sampleName/filtering/data/filtered_longreads.fa",
			"$outdir/$sampleName/$suffix/preassembly/aln.cmp.h5"
		);
		if(!$rO_jobSamtoh5->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "samtoh5", $stepName, 'SAMTOH5', $rO_jobBlasr->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobSamtoh5); 
		}

		# sort .cmp.h5
		my $rO_jobSorth5 = SmrtAnalysis::sorth5(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/preassembly/seeds.sam",
			"$outdir/$sampleName/filtering/data/filtered_longreads.fa",
			"$outdir/$sampleName/$suffix/preassembly/aln.cmp.h5"
		);
		if(!$rO_jobSorth5->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "samtoh5", $stepName, 'SAMTOH5', $rO_jobSamtoh5->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobSorth5); 
		}
		
		# quiver
		my $rO_jobQuiver = SmrtAnalysis::quiver(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/preassembly/aln.cmp.h5",
			"$outdir/$sampleName/filtering/data/filtered_longreads.fa",
			"$outdir/$sampleName/$suffix/preassembly/corrected.fa",
			"$outdir/$sampleName/$suffix/preassembly/corrected_variants.gff"
		);
		if(!$rO_jobQuiver->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "quiver", $stepName, 'QUIVER', $rO_jobSamtoh5->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobQuiver); 
		}
	 
		return $rO_jobQuiver->getCommandJobId(0);	

	}elsif($hgapAlgorithm == 2 ){
		
		## replaced pacbio smrtanalysis by all steps separetely for hgap2.
		# Separte reads using cutoff.?
		# blasr filtered_subreads.fa filtered_longreads -out seeds.m4 -m 4 -nproc 8 -bestn 24 -nCandidates 24 -noSplitSubreads -minReadLength 200 -maxScore -1000 -maxLCPLengt 16
		# echo /path/to/seeds.m4 > /path/to/seeds.m4.fofn
	
		# m4topre.py seeds.m4 /path/to/seeds.m4.fofn filtered_subreads.fa 24 > aln.pre
		# pbdagcon -a -j 8 aln.pre > corrected.fa
		# cat corrected.fa | awk 'blabla' > corrected.fq
	
		# Separate Long reads from other reads.
		my $rO_jobSplitReads = PacBioTools::splitReads(
			$rH_cfg,
			"$outdir/$sampleName/filtering/data/filtered_subreads.fasta",
			"$outdir/$sampleName/$suffix/preassemblyMinReadSize.txt",
			"$outdir/$sampleName/filtering/data/filtered_shortreads.fa",	
			"$outdir/$sampleName/filtering/data/filtered_longreads.fa"
		);
		if(!$rO_jobSplitReads->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "splitReads", $stepName, 'SPLITREADS', $dependency, $sampleName."_".$suffix, $rO_jobSplitReads); 
		}	
	
		# blasr
		my $rO_jobBlasr = SmrtAnalysis::blasr(
			$rH_cfg,
			"$outdir/$sampleName/filtering/data/filtered_subreads.fasta",
			"$outdir/$sampleName/filtering/data/filtered_longreads.fa",
			"$outdir/$sampleName/$suffix/preassembly/seeds.m4",
			"$outdir/$sampleName/$suffix/preassembly/seeds.m4.fofn",
			"nosam"
		);
		if(!$rO_jobBlasr->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "blasr", $stepName, 'BLASR', $rO_jobSplitReads->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobBlasr); 
		}
		
		# m4topre
		my $rO_jobM4topre = SmrtAnalysis::m4topre(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/preassembly/seeds.m4",
			"$outdir/$sampleName/$suffix/preassembly/seeds.m4.fofn",
			"$outdir/$sampleName/filtering/data/filtered_subreads.fasta",
			"$outdir/$sampleName/$suffix/preassembly/aln.pre"
		);
		if(!$rO_jobM4topre->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "m4topre", $stepName, 'M4TOPRE', $rO_jobBlasr->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobM4topre); 
		}
		
		# pbdagcon
		my $rO_jobPbdagcon = SmrtAnalysis::pbdagcon(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/preassembly/aln.pre",
			"$outdir/$sampleName/$suffix/preassembly/corrected.fa",
			"$outdir/$sampleName/$suffix/preassembly/corrected.fq"
		);
		if(!$rO_jobPbdagcon->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "pbdagcon", $stepName, 'PBDAGCON', $rO_jobM4topre->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobPbdagcon); 
		}
	 
		return $rO_jobPbdagcon->getCommandJobId(0);	
	}
}

#################
# Assembly (CELERA)
#
sub assembly{
	
	my $stepName	= shift;
	my $sampleName 	= shift;
	my $suffix	 	= shift; #will be ?X coverage value.
	my $merSize     = shift;
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;
	
	my $cmd;

	my $merSizeValue = $merSize;
	$merSize = "merSize".$merSize;
	
	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize/assembly\n";
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/assembly") ;#if$noMsub;

	# TODO Kmer Genie?

	# Celera config.
	my $rO_jobCeleraConfig = PacBioTools::celeraConfig(
		$rH_cfg,
		$merSizeValue,
		LoadConfig::getParam($rH_cfg, 'default', 'celeraSettings'),
		"$outdir/$sampleName/$suffix/$merSize/celeraAssembly.ini"
	);
	if(!$rO_jobCeleraConfig->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "celeraConfig", $stepName , 'CELERACONFIG', $dependency, $sampleName."_".$suffix."_$merSize", $rO_jobCeleraConfig); 
	}

	# Convert fastq to celera
	my $rO_jobFastqToCA = Celera::fastqToCA(
		$rH_cfg,
		$sampleName."_".$suffix."_".$merSize,
		"$outdir/$sampleName/$suffix/preassembly/data/corrected.fastq",
		"$outdir/$sampleName/$suffix/preassembly/data/corrected.frg"
	);
	if(!$rO_jobFastqToCA->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "fastqToCA", $stepName , 'FASTQTOCA', $rO_jobCeleraConfig->getCommandJobId(0),  $sampleName."_".$suffix."_$merSize", $rO_jobFastqToCA); 
	}

	# Run Celera	
	my $rO_jobCelera = Celera::runCelera(
		$rH_cfg,
		"$outdir/$sampleName/$suffix/$merSize/assembly/",
		$sampleName."_".$suffix."_".$merSize,
		"$outdir/$sampleName/$suffix/$merSize/celeraAssembly.ini",
		"$outdir/$sampleName/$suffix/preassembly/data/corrected.frg"
	);
	if(!$rO_jobCelera->isUp2Date()){
		my $rO_jobCelera = SubmitToCluster::printSubmitCmd($rH_cfg, "celeraAssembly", $stepName , 'ASSEMBLY', $rO_jobFastqToCA->getCommandJobId(0),  $sampleName."_".$suffix."_$merSize", $rO_jobCelera); 
	}
		
	return $rO_jobCelera->getCommandJobId(0);	
}

#################
# Polishing
# BLASR followed by Quiver. Done using smrtpipe.
#
sub polishing{
	
	my $stepName	= shift;
	my $sampleName 	= shift;
	my $suffix	 	= shift; #will be ?X coverage value.
	my $merSize     = shift;
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;

	my $merSizeValue = $merSize;
	$merSize = "merSize".$merSize;
	
	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize/polishing\n";
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/polishing") ;#if$noMsub;
	
	my $cmd;
	
	# create PacBio input.xml
	my $rO_jobFofn = SmrtAnalysis::fofns(
		$rH_cfg,
		"$outdir/fofns/$sampleName.fofn",
		"$outdir/$sampleName/$suffix/$merSize/polishing/input.xml",
		"$outdir/$sampleName/$suffix/$merSize/polishing/input.fofn"
	);
	if(!$rO_jobFofn->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "fofn", $stepName , 'FOFNS', $dependency, $sampleName."_".$suffix."_$merSize", $rO_jobFofn); 
	}

	# Upload reference.
	my $rO_jobRefUpload = SmrtAnalysis::referenceUploader(
		$rH_cfg,
		"$outdir/$sampleName/$suffix/$merSize/assembly/",
		$sampleName.$suffix.$merSize,
		"$outdir/$sampleName/$suffix/$merSize/assembly/9-terminator/H41E5.ctg.fasta"	
	);
	if(!$rO_jobRefUpload->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "referenceUpload", $stepName , 'REFUPLOAD', $rO_jobFofn->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobRefUpload); 
	}

	# Prepare xml
	$cmd = "cat ".LoadConfig::getParam($rH_cfg, 'default', 'polishingSettings');
	$cmd .= " | sed \'s|REFERENCE|".$outdir."/".$sampleName."/".$suffix."/".$merSize."/assembly/".$sampleName.$suffix.$merSize." | g\' > ".$outdir."/".$sampleName."/".$suffix."/".$merSize."/polishing/polishing.xml";
	my $rO_jobRunCommand = SmrtAnalysis::runCommand(
		$rH_cfg,
		$cmd
	);
	if(!$rO_jobRunCommand->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "XML", $stepName , 'POLISHXML', $rO_jobRefUpload->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobRunCommand); 
	}

	# Run Quiver and motif detection
	my $rO_jobPolishing = SmrtAnalysis::run(
		$rH_cfg,
		"$outdir/$sampleName/$suffix/$merSize/polishing/polishing.xml",
		#"$outdir/$sampleName/filtering/input.xml",
		"$outdir/$sampleName/$suffix/preassembly/input.xml",
		"$outdir/$sampleName/$suffix/$merSize/polishing",
		"$outdir/$sampleName/$suffix/$merSize/polishing/smrtpipe.log"
	);
	if(!$rO_jobPolishing->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "smrtanalysis", $stepName , 'PREASSEMBLY', $rO_jobRunCommand->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobPolishing); 
	}
	return $rO_jobPolishing->getCommandJobId(0);	
	
##########
#compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=8  --noXML --h5mode=w --h5fn=/lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 --seed=1 --minAccuracy=0.75 --minLength=50 --useQuality  -x -minMatch 12 -x -bestn 10 -x -minPctIdentity 70.0 --placeRepeatsRandomly --tmpDir=/lb/scratch/jtrembla/tmp --debug --regionTable=/lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/post_control_regions.fofn "/lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/input.fofn" "/lb/scratch/jtrembla/pacbio_test/H41E5/101X/assembly/H41E5101X" || exit $?; echo "Task 0 completed at `date -u`" || exit $?;echo 'Alignment Complete' || exit $?; echo "Task 1 completed at `date -u`" || exit $?;loadPulses /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/input.fofn /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread || exit $?; echo "Task 2 completed at `date -u`" || exit $?;echo 'LoadPulses Complete' || exit $?; echo "Task 3 completed at `date -u`" || exit $?;
#summarizeCoverage.py --reference /lb/scratch/jtrembla/pacbio_test/H41E5/101X/assembly/H41E5101X --numRegions=500 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/alignment_summary.gff || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#gffToBed.py --name=meanCoverage --description="Mean coverage of genome in fixed interval regions" coverage /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/alignment_summary.gff > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/coverage.bed || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#loadSequencingChemistryIntoCmpH5.py --xml /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/chemistry_mapping.xml --h5 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 || exit $?; echo "Task 0 completed at `date -u`" || exit $?;echo 'Chemistry Load Complete' || exit $?; echo "Task 1 completed at `date -u`" || exit $?;
#((which h5repack && (h5repack -f GZIP=1 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5_TMP && mv /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5_TMP /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5)) || echo 'no h5repack found, continuing w/out') || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#pbsamtools.py --bam --outfile /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.sam --refrepos /lb/scratch/jtrembla/pacbio_test/H41E5/101X/assembly/H41E5101X --readGroup movie /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#cmph5tools.py -d sort --deep --inPlace /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 || exit $?; echo "Task 0 completed at `date -u`" || exit $?;echo 'Sorting Complete' || exit $?; echo "Task 1 completed at `date -u`" || exit $?;
#extractUnmappedSubreads.py /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/filtered_subreads.fasta /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/control_reads.cmp.h5 > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/unmappedSubreads.fasta || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#variantCaller.py -P/sb/programs/analyste/software/smrtanalysis/analysis/etc/algorithm_parameters/2013-05 -v -j8 --algorithm=quiver /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 -r /lb/scratch/jtrembla/pacbio_test/H41E5/101X/assembly/H41E5101X/sequence/H41E5101X.fasta -o /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff -o /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/consensus.fasta.gz -o /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/consensus.fastq.gz || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#ummarizeConsensus.py --variantsGff /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/alignment_summary.gff --output /lb/scratch/jtrembla/tmp/tmppIn9qP.gff || exit $?; echo "Task 0 completed at `date -u`" || exit $?;mv /lb/scratch/jtrembla/tmp/tmppIn9qP.gff /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/alignment_summary.gff || exit $?; echo "Task 1 completed at `date -u`" || exit $?;
#gffToBed.py --name=variants --description='PacBio: snps, insertions, and deletions derived from consensus calls against reference' variants /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.bed || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#gffToVcf.py --globalReference=H41E5101X /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.vcf || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#gzip -f /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
#########		
}

#################
# BLAST
#
sub blast{

	my $stepName	= shift;
	my $sampleName 	= shift;
	my $suffix	 	= shift; #will be ?X coverage value.
	my $merSize     = shift;
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;

	my $merSizeValue = $merSize;
	$merSize = "merSize".$merSize;
	
	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize/blast\n";
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/blast") ;#if$noMsub;

	my $cmd;
	
	# Blast contigs against nt
	my $rO_jobBlast = BLAST::dcmegablast(
		$rH_cfg,
		"gunzip -c $outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta.gz > $outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta",
		"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta",
		"7",
		"$outdir/$sampleName/$suffix/$merSize/blast/blast_report.csv"
	);
	if(!$rO_jobBlast->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "blast", "dc-megablast" , 'BLAST', $dependency,  $sampleName."_".$suffix."_$merSize", $rO_jobBlast); 
	}

	# Get fasta file of best hit.
	my $rO_jobBlastDb = BLAST::blastdbcmd(
		$rH_cfg,
		"\$(head -n 6 $outdir/$sampleName/$suffix/$merSize/blast/blast_report.csv | tail -n 1 |  awk -F \\\\t '{print \$2}' | sed 's/gi|\\([0-9]*\\)|.*/\\1/' | tr '\\n' '  ')",
		"$outdir/$sampleName/$suffix/$merSize/blast/nt_reference.fasta"
	);
	if(!$rO_jobBlastDb->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "blastdbcmd", "dc-megablast" , 'BLASTDBCMD', $rO_jobBlast->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobBlastDb); 
	}
	
	return $rO_jobBlastDb->getCommandJobId(0);	
}

#################
# MUMMER
#
sub mummer{
	
	my $stepName	= shift;
	my $sampleName 	= shift;
	my $suffix	 	= shift; #will be ?X coverage value.
	my $merSize     = shift;
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;

	my $merSizeValue = $merSize;
	$merSize = "merSize".$merSize;

	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize/mummer\n";
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/mummer") ;#if$noMsub;

	my $cmd;

	# Run nucmer
	my $rO_jobNucmer = Mummer::nucmer(
		$rH_cfg,
		100,
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer",
		"$outdir/$sampleName/$suffix/$merSize/blast/nt_reference.fasta",
		"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta"	
	);
	if(!$rO_jobNucmer->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "nucmer", $stepName , 'NUCMER', $dependency, $sampleName."_".$suffix."_$merSize", $rO_jobNucmer); 
	}

	# Generate plot.
	my $rO_jobMummerPlot = Mummer::mummerPlot(
		$rH_cfg,
		"$sampleName/$suffix-nucmer/",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.delta",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.delta"
	);
	if(!$rO_jobMummerPlot->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "mummer-plot", $stepName , 'MUMMER_PLOT', $rO_jobNucmer->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobMummerPlot); 
	}

	# Dnadiff
	my $rO_jobDnadiff = Mummer::dnadiff(
		$rH_cfg,
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.dnadiff",
		"$outdir/$sampleName/$suffix/$merSize/blast/nt_reference.fasta",
		"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta"
	);
	if(!$rO_jobDnadiff->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "dnadiff", $stepName , 'DNADIFF', $rO_jobMummerPlot->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobDnadiff); 
	}
	
	# showsnp
	my $rO_jobShowsnp = Mummer::showsnp(
		$rH_cfg,
		200,
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.dnadiff.delta",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.dnadiff.delta.snpflank"		
	);
	if(!$rO_jobShowsnp->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "showsnp", $stepName , 'SHOWSNP', $rO_jobDnadiff->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobShowsnp); 
	}

	# Run on self	
	my $rO_jobNucmerSelf = Mummer::nucmer(
		$rH_cfg,
		100,
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.self",
		"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta",
		"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta"
	);
	if(!$rO_jobNucmerSelf->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "nucmer", $stepName , 'NUCMER_SELF', $rO_jobShowsnp->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobNucmerSelf); 
	}
	
	# Make plot of self comparison.
	my $rO_jobMummerPlotSelf = Mummer::mummerPlot(
		$rH_cfg,
		"$sampleName/$suffix-nucmer-self/",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.self.delta",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.self.delta"
	);
	if(!$rO_jobMummerPlotSelf->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "mummer-plot-self", $stepName, 'MUMMER_PLOT_SELF', $rO_jobNucmerSelf->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobMummerPlotSelf); 
	}

	return $rO_jobMummerPlotSelf->getCommandJobId(0);
}

#################
# Epigenome TODO
#
sub epigenome{
	
	my $stepName	= shift;
	my $sampleName 	= shift;
	my $suffix	 	= shift; #will be ?X coverage value.
	my $merSize     = shift;
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;

	my $merSizeValue = $merSize;
	$merSize = "merSize".$merSize;

	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize/motifs\n";
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/motifs") ;#if$noMsub;

	my $cmd;

	# Prepare xml
	$cmd = "cat ".LoadConfig::getParam($rH_cfg, 'default', 'motifsSettings');
	$cmd .= " | sed \'s|REFERENCE|".$outdir."/".$sampleName.$suffix."/assembly/".$sampleName.$suffix." | g\' > ".$outdir."/".$sampleName.$suffix."/motifs/motifs.xml";
	my $rO_jobRunCommand = SmrtAnalysis::runCommand(
		$rH_cfg,
		$cmd
	);
	if(!$rO_jobRunCommand->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "XML", $stepName , 'MOTIFSXML', $dependency, $sampleName."_".$suffix, $rO_jobRunCommand); 
	}

	# smart pipe to create fastq in --output <dir> arg.
	my $rO_jobMotifs = SmrtAnalysis::run(
		$rH_cfg,
		"$outdir/$sampleName/$suffix/$merSize/motifs.xml",
		"$outdir/$sampleName/$suffix/$merSize/motifs/input.xml",
		"$outdir/$sampleName/$suffix/$merSize/motifs",
		"$outdir/$sampleName/$suffix/$merSize/motifs/smrtpipe.log"
	);
	if(!$rO_jobMotifs->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "smrtanalysis", $stepName , 'MOTIFS', $rO_jobRunCommand->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobMotifs); 
	}
	
	return $rO_jobMotifs->getCommandJobId(0); 
}

#################
# Report
#
sub report {

	my $stepName	= shift;
	my $sampleName 	= shift;
	my $suffix	 	= shift; #will be ?X coverage value.
	my $merSize     = shift;
	my $filePath	= shift;
	my $rH_cfg 		= shift;
	my $dependency 	= shift;

	my $merSizeValue = $merSize;
	$merSize = "merSize".$merSize;

	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize/report\n";
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize") ;#if$noMsub;
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/report") ;#if$noMsub;

	my $cmd;

	# Generate table(s) and figures.
	my $rO_jobAssemblyStats = PacBioTools::assemblyStats(
		$rH_cfg,
		#"$outdir/$sampleName/$suffix/",
		"$outdir/$sampleName/$suffix/preassembly/data/filtered_summary.csv",
		"$outdir/$sampleName/$suffix/preassembly/results/filterReports_filterStats.xml",
		"$outdir/$sampleName/$suffix/$merSize/assembly/9-terminator/$sampleName.qc",
		$sampleName,
		$suffix."_".$merSize,
		$estimatedGenomeSize,
		$smrtCells,
		"$outdir/$sampleName/$suffix/$merSize/report"
	);
	if(!$rO_jobAssemblyStats->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "assembly-stats", $stepName, 'ASSEMBLY_STATS', $dependency, $sampleName, $rO_jobAssemblyStats); 
	}

	# Generate report with Noozle
	my $title = ""; 
	my $titleTMP = LoadConfig::getParam($rH_cfg, 'report','projectName');
	if (defined($titleTMP) && !($titleTMP eq "")) {
		$title= 'report.title=\"' .$titleTMP .'\",';
	}   
	#my $path = ""; 
	#my $pathTMP = LoadConfig::getParam($rH_cfg, 'report','report.path');
    #if (defined($pathTMP) && !($pathTMP eq "")) {
	#	$path= 'report.path=\"' .$pathTMP .'\",';
    #}   
	my $author = ""; 
	my $authorTMP = LoadConfig::getParam($rH_cfg, 'report','report.author');
	if (defined($authorTMP) && !($authorTMP eq "")) {
		$author= 'report.author=\"' .$authorTMP .'\",';
	}   
	my $contact = ""; 
	my $contactTMP = LoadConfig::getParam($rH_cfg, 'report','report.contact');
	if (defined($contactTMP) && !($contactTMP eq "")) {
		$contact= 'report.contact=\"' .$contactTMP .'\",';
	} 	
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs(undef, undef);

	$cmd .= 'module load ' .LoadConfig::getParam($rH_cfg, 'report','moduleVersion.cranR') .' &&';
	$cmd .= ' R --vanilla -e \'library(gqSeqUtils) ;';
	$cmd .= ' mugqicPipelineReport(';
	$cmd .= ' pipeline=\"PacBioAssembly\",';
	$cmd .= ' ' .$title;
	$cmd .= ' ' .$author;
	$cmd .= ' ' .$contact;
	$cmd .= ' ini.file.path=\"' . $configFile . '\",';
	$cmd .= ' project.path=\"' . "$outdir/$sampleName/$suffix/$merSize" . '\",';
	$cmd .= ' report.path=\"' . "$outdir/$sampleName/$suffix/$merSize/report" . '\")\'';

	# Here I had could not use the GqSeqUtils class because it does not support a project.path arg...	
	if (!$ro_job->isUp2Date()) {
		$ro_job->addCommand($cmd);
	}

	if (!$ro_job->isUp2Date()) {

		# Build command string to submit to cluster.
	    SubmitToCluster::printSubmitCmd(
	        $rH_cfg, 
			$stepName,
			$sampleName,
			$suffix,
			$rO_jobAssemblyStats->getCommandJobId(0),
			$sampleName,
			$ro_job
		);
	}	

	return $ro_job->getCommandJobId(0);
}

1;
