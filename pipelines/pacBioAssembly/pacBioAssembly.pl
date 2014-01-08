#!/usr/bin/env perl

=head1 NAME

I<pacBioAssembly>

=head1 SYNOPSIS

pacBioAssembly.pl

=head1 DESCRIPTION

B<pacBioAssembly> Is the main PacBio HGAP assembly pipeline.

=head1 AUTHOR

B<Julien Tremblay> - I<julien.tremblay@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing

=head1 SUBROUTINES

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
my $smrtCells;
my $samplesFile;
my $noMsub;
my $estimatedGenomeSize;
my $estimatedCoverage;
my $numberOfBases;
my $libType;
my $verbose = 0;
my $TMPDIR;
my $hgapAlgorithm;
my $computeCutoffFlag;

&main();

sub printUsage {
  print "\nUsage: perl ".$0." \n";
  print "\t-c  config file\n";
  print "\t-s  start step, inclusive\n";
  print "\t-e  end step, inclusive\n";
  print "\t-n  nanuq sample sheet\n";
  print "\t-m  flag to run job interactively\n";
  print "\t-a  0: hgap2 legacy(smrtpipe). 1: hgap2, separated commands (Default)\n";
  print "\t-t  0: Do not compute minReadLength threshold for preassembly (i.e. use smrtpipe default). 1: Compute minReadLength threshold for preassembly. Default=1. Do not use only for development.\n";
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
	getopts('c:s:e:n:w:m:v:a:t:', \%opts);
	
	if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'e'}) || !defined($opts{'n'}) ) {
		printUsage();
		exit(1);
	}
	
	# Load config file
	$noMsub = 1 if(defined($opts{'m'}));
	if(!defined($opts{'t'})){
		$computeCutoffFlag = 1 
	}elsif($opts{'t'} == 0){
		$computeCutoffFlag = 0;
	}elsif($opts{'t'} == 1){
		$computeCutoffFlag = 1;
	}

	my %jobIdVarPrefix;
	my %cfg 				= LoadConfig->readConfigFile($opts{'c'});
	$TMPDIR				    = LoadConfig::getParam(\%cfg, 'default', 'tmpDir');	

	# Get design groups
	$configFile 			= abs_path($opts{'c'});
	$samplesFile			= abs_path($opts{'n'}); 
	$hgapAlgorithm          = $opts{'a'};
	$hgapAlgorithm          = 1 if(!defined($hgapAlgorithm));
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
	my %hSeen;

	open(SAMPLES, '<'.$samplesFile) or die "Can't find file $samplesFile\n";
	while(<SAMPLES>){
		chomp;
		my @row = split(/\t/, $_);
		$hohSamples{$row[0]}{$row[1]}{$row[2]} = $row[2];
		$hLibType{$row[0]} = $row[3];
			if($row[2] =~ m/\.bas\.h5/){
				$hNumberOfBases{$row[0]} = $hNumberOfBases{$row[0]} + $row[4];
			}else{
				$row[2] =~ m/\/(.*)\.\d+\.bax\.h5/;
				$hSeen{$1}++;
				if(exists $hNumberOfBases{$row[0]}){
					$hNumberOfBases{$row[0]} = $hNumberOfBases{$row[0]} + $row[4] if($hSeen{$1} == 1);
				}else{
					$hNumberOfBases{$row[0]} = $row[4] if($hSeen{$1} == 1);
				}
				$hSmrtCells{$row[0]}++ if($hSeen{$1} == 1);
			}
		$hEstimatedGenomeSize{$row[0]} = $row[5];
	}
	close(SAMPLES);	
	
	# fofns will be written here. One fofn per sample for bas5 (optional), 
	# but 3 per samples for bax. One smrt cell can the same sample more than one time.
	my %fofns;
	while(my($sampleName, $value1)=each(%hohSamples)){
		#print STDERR "$sampleName\n";		
		
		# Write fofns here in Main
		open(FOFN, ">$outdir/fofns/$sampleName.fofn") or die "Can't open fofn file for writing ... $!\n";
		
		while(my($well, $value2)=each(%$value1)){
			while(my($pathToFofn, $value3)=each(%$value2)){
				#print STDERR $sampleName."\t".$well."\t".$pathToFofn."\n";
				print FOFN "$value3\n";
			}
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
		$estimatedCoverage = int($numberOfBases / $estimatedGenomeSize);
		
		print STDERR "[DEBUG] EstimatedGenomeSize: ".$estimatedGenomeSize."\n";
		print STDERR "[DEBUG] smrtCells: ".$smrtCells."\n";
		print STDERR "[DEBUG] LibType: ".$libType."\n";
		print STDERR "[DEBUG] numberOfbases: ".$numberOfBases."\n";
		print STDERR "[DEBUG] Estimated coverage: ".$estimatedCoverage."\n";	

		my $currentStep;
		my $rootDependency = undef;
		my $dependency = undef;

		# global.
		for($currentStep = $opts{'s'}-1; $currentStep <= ($opts{'e'}-1); $currentStep++) {

  			my $fname = $steps[$currentStep]->{'name'};
   			my $subref = \&$fname;
     
			if ($steps[$currentStep]->{'stepLoop'} eq 'global') {
				# Tests for the first step in the list. Used for dependencies.
				$dependency = &$subref($steps[$currentStep]->{'name'}, $currSampleName, $currFofn, \%cfg, $rootDependency); # Make sure only one sub/step is labeled as 'global' 
				$globalDep{$fname}->{$currSampleName} = $dependency;
			}   
		}   
		
		# do assembly for the various coverage values entered in the ini file. And also do assembly with
		# the "logical" total sequenced bp / genome estimated size.	
		my $coverageRange = LoadConfig::getParam(\%cfg, 'default', 'coverageCutoff'); 
		my @coverageRange = split(/:/, $coverageRange);
		#push(@coverageRange, $estimatedCoverage);
		
		my $merSizes = LoadConfig::getParam(\%cfg, 'default', 'merSizes');
		my @merSizes = split(/:/, $merSizes);

		# Sample/assembly
		foreach my $coverage(@coverageRange){
		
			my $dependency2 = $dependency;
	
			for($currentStep = $opts{'s'}-1; $currentStep <= ($opts{'e'}-1); $currentStep++) {
	
  				my $fname = $steps[$currentStep]->{'name'};
   				my $subref = \&$fname;
   	  
				if ($steps[$currentStep]->{'stepLoop'} eq 'sample') { #preassembly
					if($steps[$currentStep]->{'name'} eq 'getStats'){
						$dependency2 = &$subref($steps[$currentStep]->{'name'}, $currSampleName, $estimatedCoverage, $coverage, $currFofn, \%cfg, $dependency2);
						$globalDep{$fname}->{$currSampleName} = $dependency;

					}else{
						# Tests for the first step in the list. Used for dependencies.
						$dependency2 = &$subref($steps[$currentStep]->{'name'}, $currSampleName, $coverage."percent", $currFofn, \%cfg, $dependency2);
						$globalDep{$fname}->{$currSampleName} = $dependency;
					}
				}elsif($steps[$currentStep]->{'stepLoop'} eq 'assembly') {
			
					my $dependency3 = $dependency2;
				
					foreach my $merSize(@merSizes){			
						print STDERR "[DEBUG] MERSIZE: ".$merSize."\n";
						# Tests for the first step in the list. Used for dependencies.
						$dependency2 = &$subref($steps[$currentStep]->{'name'}, $currSampleName, $coverage."percent", $merSize, $currFofn, \%cfg, $dependency2); 
						$globalDep{$fname}->{$currSampleName} = $dependency;
					}
				} 
			}  	
		}
	}
	exit;		
}

=head2 filtering()

Filtering. This step uses smrtpipe.py (From the SmrtAnalysis package) and will filter reads and subreads  based on their length and QVs.
1- fofnToSmrtpipeInput.py. 
2- modify RS_Filtering.xml files according to reads filtering values entered in .ini file.
3- smrtpipe.py with filtering protocol
4- prinseq-lite.pl to write fasta file based on fastq file.
Informative run metrics such as loading efficiency, readlengths, and base quality are generated in this step as well..

=cut

############
# Filtering
#
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
	system("mkdir -p $outdir/$sampleName");
	system("mkdir -p $outdir/$sampleName/filtering");
	
	# Define tmpdir.
	my $tmpdir = LoadConfig::getParam($rH_cfg, 'smrtanalysis', 'tmpDir')."/".$sampleName."_".$stepName;
	print STDOUT "mkdir -p ".$tmpdir."\n";

	# smart pipe to create fastq in --output <dir> arg.
	my $rO_jobFiltering = SmrtAnalysis::filtering(
		$rH_cfg,
		# Fofn
		"$outdir/fofns/$sampleName.fofn", 
		"$outdir/$sampleName/filtering/input.xml",
		"$outdir/$sampleName/filtering/input.fofn",	
		# Xml
		LoadConfig::getParam($rH_cfg, 'default', 'filteringSettings'),
		$outdir."/".$sampleName."/filtering.xml",		
		# Filtering smrtpipe
		"$outdir/$sampleName/filtering",
		"$outdir/$sampleName/filtering/smrtpipe.log",
		$tmpdir
	);
	if(!$rO_jobFiltering->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "filtering", $stepName , 'FILTERING_SMRTPIPE', $dependency, $sampleName, $rO_jobFiltering); 
	}

 	return $rO_jobFiltering->getCommandJobId(0);	
}

=head2 getStats()

Cutoff value for splitting long reads from short reads is done here using 
estimated coverage and estimated genome size.

You should estimate the overall coverage and length distribution for putting in
the correct options in the configuration file. You will need to decide a
length cutoff for the seeding reads. The optimum cutoff length will depend on
the distribution of the sequencing read lengths, the genome size and the
overall yield. Here, you provide a percentage value that corresponds to the 
fraction of coverage you want to use as seeding reads.

First, loop through fasta sequences. put the length of each sequence in an array, 
sort it, loop through it again and compute the cummulative length coveredby each 
sequence as we loop though the array. Once that length is > (coverage * genome 
size) * $percentageCutoff (e.g. 0.10), we have our threshold. The idea is to 
consider all reads above that threshold to be seeding reads to which will be 
align lower shorter subreads.

=cut

#################
# getStats
#
sub getStats{
	
	my $stepName	      = shift;
	my $sampleName 	      = shift;
	my $estimatedCoverage = shift; #will be ?X coverage value.
	my $coverageCutoff    = shift;
	my $filePath	      = shift;
	my $rH_cfg 		      = shift;
	my $dependency 	      = shift;

	my $suffix = $coverageCutoff."percent";

	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/preassembly\n";
	system("mkdir -p $outdir/$sampleName/$suffix");
	system("mkdir -p $outdir/$sampleName/$suffix/preassembly");

	$coverageCutoff = $coverageCutoff / 100; 
	
	# Choose a subread length threshold such that subreads above the threshold provide about 20x coverage of the genome.
	my $rO_jobGetCutoff = PacBioTools::getCutoff(
		$rH_cfg,
		"$outdir/$sampleName/filtering/data/filtered_subreads.fasta",
		$estimatedCoverage,
		$estimatedGenomeSize,
		$coverageCutoff,
		LoadConfig::getParam($rH_cfg, 'default', 'preassemblySettings'),
		"$outdir/$sampleName/$suffix/preassembly.xml",
		"$outdir/$sampleName/$suffix/preassemblyMinReadSize.txt"
	);
	if(!$rO_jobGetCutoff->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "pacbio_tools", $stepName , 'GETCUTOFF', $dependency, $sampleName."_".$suffix, $rO_jobGetCutoff); 
	}
 
	return $rO_jobGetCutoff->getCommandJobId(0);
}

=head2 preAssembly()

Having in hand a cutoff value, filtered reads are splitted between short and long reads. Short reads
are aligned against long reads and consensus (e.g. corrected reads) are generated from these alignments.
1- split reads between long and short.
2- blasr (Aligner for PacBio reads)
3- m4topre (Converts .m4 blasr output in .pre format.)
4- pbdagcon (generates corrected reads from alignments)

=cut

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

	# Define outdir
	my $outdir = LoadConfig::getParam($rH_cfg, 'default', 'outdir');
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/preassembly\n";
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/preassembly/data\n";
	system("mkdir -p $outdir/$sampleName/$suffix");
	system("mkdir -p $outdir/$sampleName/$suffix/preassembly");
	system("mkdir -p $outdir/$sampleName/$suffix/preassembly/data");

	# Define tmpdir.
	my $tmpdir = LoadConfig::getParam($rH_cfg, 'smrtanalysis', 'tmpDir')."/".$sampleName."_".$stepName."_".$suffix;
	print STDOUT "mkdir -p ".$tmpdir."\n";

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
		my $currentParams;
		if( $computeCutoffFlag == 1 ){ # do not use default, compute cutoff.
			$currentParams = "$outdir/$sampleName/$suffix/preassembly.xml"
		}elsif( $computeCutoffFlag == 0 ){ # use default, do not compute cutoff
			$currentParams = LoadConfig::getParam($rH_cfg, 'default', 'preassemblySettings')
		}else{ die "Invalid value for -t arg...\n";}

		my $rO_jobPreassembly = SmrtAnalysis::run(
			$rH_cfg,
			$currentParams,
			"$outdir/$sampleName/$suffix/preassembly/input.xml",
			"$outdir/$sampleName/$suffix/preassembly",
			"$outdir/$sampleName/$suffix/preassembly/smrtpipe.log",
			$tmpdir
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
			"$outdir/$sampleName/filtering/data/filtered_subreads.fasta",
			"$outdir/$sampleName/filtering/data/filtered_longreads.fa",
			"$outdir/$sampleName/$suffix/preassembly/data/seeds.m4",
			"$outdir/$sampleName/$suffix/preassembly/data/seeds.m4.fofn",
			"nosam"
		);
		if(!$rO_jobBlasr->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "blasr", $stepName, 'BLASR', $rO_jobSplitReads->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobBlasr); 
		}
		
		# m4topre
		my $rO_jobM4topre = SmrtAnalysis::m4topre(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/preassembly/data/seeds.m4",
			"$outdir/$sampleName/$suffix/preassembly/data/seeds.m4.fofn",
			"$outdir/$sampleName/filtering/data/filtered_subreads.fasta",
			"$outdir/$sampleName/$suffix/preassembly/data/aln.pre"
		);
		if(!$rO_jobM4topre->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "m4topre", $stepName, 'M4TOPRE', $rO_jobBlasr->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobM4topre); 
		}
		
		# pbdagcon
		my $rO_jobPbdagcon = SmrtAnalysis::pbdagcon(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/preassembly/data/aln.pre",
			"$outdir/$sampleName/$suffix/preassembly/data/corrected.fasta",
			"$outdir/$sampleName/$suffix/preassembly/data/corrected.fastq"
		);
		if(!$rO_jobPbdagcon->isUp2Date()){
			SubmitToCluster::printSubmitCmd($rH_cfg, "pbdagcon", $stepName, 'PBDAGCON', $rO_jobM4topre->getCommandJobId(0), $sampleName."_".$suffix, $rO_jobPbdagcon); 
		}
	 
		return $rO_jobPbdagcon->getCommandJobId(0);	
	}
}

=head2 assembly()

Corrected reads are assembled to generates contigs. Please see Celera documentation.
http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=RunCA#ovlThreads
Quality of assembly seems to be highly sensitive to paramters you give Celera.
1- Generate celera config files using paramters provided in the .ini file.
2- fastqToCA. Generates input file compatible with the Celera assembler
3- runCA. Run the Celera assembler.

=cut

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
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize");
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/assembly");

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

=head2 polishing()

Align raw reads on the Celera assembly with BLASR. Load pulse information from bax files into aligned file. Sort that file
and run quiver (variantCaller.py).

1- Generate fofn
2- Upload Celera assembly with smrtpipe refUploader
3- Compare sequences
4- Load pulses
5- Sort .cmp.h5 file
6- variantCaller.py  

=cut

#################
# Polishing (BLASR + Quiver).
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
	print STDOUT "mkdir -p $outdir/$sampleName/$suffix/$merSize/polishing/data\n";
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize");
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/polishing");
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/polishing/data");
	
	my $tmpdir = LoadConfig::getParam($rH_cfg, 'smrtanalysis', 'tmpDir')."/".$sampleName."_".$stepName."_".$suffix."_".$merSize;
	print STDOUT "mkdir -p ".$tmpdir."\n";
	
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
		"$outdir/$sampleName/$suffix/$merSize/assembly/9-terminator/$sampleName"."_".$suffix."_".$merSize.".ctg.fasta"	
	);
	if(!$rO_jobRefUpload->isUp2Date()){
		SubmitToCluster::printSubmitCmd($rH_cfg, "referenceUpload", $stepName , 'REFUPLOAD', $rO_jobFofn->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobRefUpload); 
	}
	
	if($hgapAlgorithm == 0){
	
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
			"$outdir/$sampleName/$suffix/$merSize/polishing/smrtpipe.log",
			$tmpdir
		);
		if(!$rO_jobPolishing->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "smrtanalysis", $stepName , 'POLISHING', $rO_jobRunCommand->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobPolishing); 
		}
		return $rO_jobPolishing->getCommandJobId(0);	
	
	}elsif($hgapAlgorithm == 1){
		my $rO_jobCompareSequences = SmrtAnalysis::compareSequences(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/aligned_reads.cmp.h5",
			"$outdir/$sampleName/filtering/data/filtered_regions.fofn", 
			"$outdir/$sampleName/$suffix/$merSize/polishing/input.fofn",
			$outdir."/".$sampleName."/".$suffix."/".$merSize."/assembly/".$sampleName.$suffix.$merSize,
			$tmpdir
		);
		if(!$rO_jobCompareSequences->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "compareSequences", $stepName , 'POLISH', $rO_jobRefUpload->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobCompareSequences); 
		}

		my $rO_jobLoadPulses = SmrtAnalysis::loadPulses(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/aligned_reads.cmp.h5",
			"$outdir/$sampleName/$suffix/$merSize/polishing/input.fofn"
		);
		if(!$rO_jobLoadPulses->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "loadPulses", $stepName , 'LOADPULSES', $rO_jobCompareSequences->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobLoadPulses); 
		}

		my $rO_jobSortH5 = SmrtAnalysis::sortH5(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/aligned_reads.cmp.h5",
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/aligned_reads.cmp.h5.sorted",
		);
		if(!$rO_jobLoadPulses->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "sortH5", $stepName , 'SORTH5', $rO_jobLoadPulses->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobSortH5); 
		}

		my $rO_jobCallVariants = SmrtAnalysis::variantCaller(
			$rH_cfg,
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/aligned_reads.cmp.h5.sorted",  #cmpH5alignedReads
			$outdir."/".$sampleName."/".$suffix."/".$merSize."/assembly/".$sampleName.$suffix.$merSize."/sequence/".$sampleName.$suffix.$merSize.".fasta",
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/variants.gff",
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta.gz",
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fastq.gz"
		);
		if(!$rO_jobCallVariants->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "callVariants", $stepName , 'CALLVARIANTS', $rO_jobSortH5->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobCallVariants); 
		}
		
		my $rO_jobSummarizePolishing = SmrtAnalysis::summarizePolishing(
			$rH_cfg,
			$sampleName."_".$suffix."_".$merSize,
			$outdir."/".$sampleName."/".$suffix."/".$merSize."/assembly/".$sampleName.$suffix.$merSize,                       # REFERENCE
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/aligned_reads.cmp.h5.sorted",                                # cmpH5alignedReads
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/alignment_summary.gff",                                      # alignment_summary.gff 
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/coverage.bed",                                               # coverage.bed
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/aligned_reads.sam",                                          # sam file
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/variants.gff",                                               # alignment_summary2.gff 
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/variants.bed",                                               # alignment_summary2.gff 
			"$outdir/$sampleName/$suffix/$merSize/polishing/data/variants.vcf"                                                # alignment_summary2.gff 
		);
		if(!$rO_jobSummarizePolishing->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "summarizePolishing", $stepName , 'SUMMARIZE_POLISHING', $rO_jobCallVariants->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobSummarizePolishing); 
		}
		
		return $rO_jobSummarizePolishing->getCommandJobId(0);	
	}	
}

=head2 blast()

Blast polished assembly against nr using dc-megablast.

=cut

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
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize");
	system("mkdir -p $outdir/$sampleName/$suffix/$merSize/blast");
	
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
		"\\\$(head -n 6 $outdir/$sampleName/$suffix/$merSize/blast/blast_report.csv | tail -n 1 |  awk -F \\\\\\t '{print \\\$2}' | sed 's/gi|\\([0-9]*\\)|.*/\\1/' | tr '\\n' '  ')",
		"$outdir/$sampleName/$suffix/$merSize/blast/nt_reference.fasta"
	);
	if(!$rO_jobBlastDb->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "blastdbcmd", "dc-megablast" , 'BLASTDBCMD', $rO_jobBlast->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobBlastDb); 
	}
	
	return $rO_jobBlastDb->getCommandJobId(0);	
}

=head2 mummer()

Using MUMmer, align polished assembly against best hit from blast job. Also align polished assembly against itself
to detect structure variation such as repeats, etc.

=cut

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

	# Run nucmer
	my $rO_jobMummerRef = Mummer::reference(
		$rH_cfg,
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer",
		"$outdir/$sampleName/$suffix/$merSize/blast/nt_reference.fasta",
		"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta",
		$sampleName."_".$suffix."-nucmer_".$merSize,
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.delta",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.delta",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.dnadiff",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.dnadiff.delta",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.dnadiff.delta.snpflank"
	);
	if(!$rO_jobMummerRef->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "mummerRef", $stepName , 'MUMMER_REF', $dependency, $sampleName."_".$suffix."_$merSize", $rO_jobMummerRef); 
	}

	my $rO_jobMummerSelf = Mummer::self(
		$rH_cfg,
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.self",
		"$outdir/$sampleName/$suffix/$merSize/polishing/data/consensus.fasta",
		$sampleName."_".$suffix."-nucmer-self_".$merSize,
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.self.delta",
		"$outdir/$sampleName/$suffix/$merSize/mummer/$sampleName.nucmer.self.delta"
	);
	if(!$rO_jobMummerSelf->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "mummerSelf", $stepName , 'MUMMER_SELF', $rO_jobMummerRef->getCommandJobId(0), $sampleName."_".$suffix."_$merSize", $rO_jobMummerSelf); 
	}
	
	return $rO_jobMummerSelf->getCommandJobId(0);
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


=head2 report()

Generates summary tables and Generates MUGQIC style nozzle report.

=cut

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
		"$outdir/$sampleName/filtering/data/filtered_summary.csv",
		"$outdir/$sampleName/$suffix/$merSize/assembly/9-terminator/".$sampleName."_".$suffix."_".$merSize.".qc",
		"$outdir/$sampleName/$suffix/$merSize/assembly/9-terminator/$sampleName"."_".$suffix."_".$merSize.".ctg.fasta",
		$sampleName,
		$suffix."_".$merSize,
		$estimatedGenomeSize,
		$smrtCells,
		"$outdir/$sampleName/$suffix/$merSize/report"
	);
	if(!$rO_jobAssemblyStats->isUp2Date()) {
		SubmitToCluster::printSubmitCmd($rH_cfg, "assembly-stats", $stepName, 'ASSEMBLY_STATS', $dependency, $sampleName."_".$suffix."_$merSize", $rO_jobAssemblyStats); 
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
	$cmd .= ' R --no-save -e \'library(gqSeqUtils) ;';
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
			"NOZZLE_".$sampleName."_".$suffix."_".$merSize,
			$rO_jobAssemblyStats->getCommandJobId(0),
			$sampleName."_".$suffix."_$merSize",
			$ro_job
		);
	}	

	return $ro_job->getCommandJobId(0);
}

1;
