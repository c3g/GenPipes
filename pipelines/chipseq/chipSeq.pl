#!/usr/bin/perl

=head1 NAME

I<chipSeq>

=head1 SYNOPSIS

chipSeq.pl

=head1 DESCRIPTION

B<chipSeq> is the main ChIPseq pipeline.

=head1 AUTHOR

B<Johanna Sandoval> - I<johanna.sandoval@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debug

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
}


# Dependencies
#--------------------
use Getopt::Std;

use LoadConfig;
use SampleSheet;
use SAMtools;
use Picard;
use HtseqCount;
use SequenceDictionaryParser;
use SubmitToCluster;
use Homer;
use MACS2;
use BWA;
use Trimmomatic;
use Metrics;
use Wiggle;

#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'trimming' , 'stepLoop' => 'sample' , 'output' => 'reads'});
push(@steps, {'name' => 'aligning' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
push(@steps, {'name' => 'rawCounts' , 'stepLoop' => 'sample' , 'output' => 'raw_counts'});
push(@steps, {'name' => 'qcTagDirectories' , 'stepLoop' => 'sample' , 'output' => 'tags'});
push(@steps, {'name' => 'wiggle' , 'stepLoop' => 'sample' , 'output' => 'tracks'});
push(@steps, {'name' => 'peakCall' , 'stepLoop' => 'design' , 'output' => 'peak_call'});
push(@steps, {'name' => 'annotation' , 'stepLoop' => 'design' , 'output' => 'annotation'});
push(@steps, {'name' => 'motif' , 'stepLoop' => 'design' , 'output' => 'motif'});
push(@steps, {'name' => 'qcPlots' , 'stepLoop' => 'general' , 'output' => 'graphs'});
# push(@steps, {'name' => 'report' , 'stepLoop' => 'group' , 'output' => 'DGE'});
# push(@steps, {'name' => 'deliverable' , 'stepLoop' => 'group' , 'output' => 'Deliverable'});


my %globalDep;
for my $stepName (@steps) { 
	$globalDep{$stepName -> {'name'} } ={};
}

my $designFilePath;
my $workDirectory;


&main();

sub printUsage {
  print "\nUsage: perl ".$0." \n";
  print "\t-c  config file\n";
  print "\t-s  start step, inclusive\n";
  print "\t-e  end step, inclusive\n";
  print "\t-n  nanuq sample sheet\n";
  print "\t-d  design file\n";
  print "\t-w  work directory\n";
  print "\n";
  print "Steps:\n";
  for(my $idx=0; $idx < @steps; $idx++) {
    print "".($idx+1).'- '.$steps[$idx]->{'name'}."\n";
  }
  print "\n";
}

sub main {
	my %opts;
	getopts('c:s:e:n:d:w:', \%opts);
	
	if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'e'}) || !defined($opts{'n'}) || !defined($opts{'d'}) || !defined($opts{'w'} ) ) {
		printUsage();
		exit(1);
	}
	
	my %jobIdVarPrefix;
	my %cfg = LoadConfig->readConfigFile($opts{'c'});
	my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($opts{'n'});
	my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);
	$designFilePath = $opts{'d'};
	# get Design groups
	my $rHoAoA_designGroup = MACS2::getDesign(\%cfg,$designFilePath);
	$workDirectory = $opts{'w'};
	#generate sample jobIdprefix
	my $cpt = 1;
	
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
		my $cpt2=1;
		$jobIdVarPrefix{$sampleName} = $cpt;
		my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
		for my $rH_laneInfo (@$rAoH_sampleLanes) {
			$jobIdVarPrefix{$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} = $cpt.'_'.$cpt2;
			$cpt2++;
		}
		$cpt++;
	}
	$cpt = 1;
	# generate design jobIdprefix
	for my $designName (keys %{$rHoAoA_designGroup}) {
		$jobIdVarPrefix{$designName} = $cpt;
		$cpt++;
	}
	
	print "cd $workDirectory \n";

	my $latestBam;
	
	for(my $current = $opts{'s'}-1; $current <= ($opts{'e'}-1); $current++) {
		my $fname = $steps[$current]->{'name'};
		my $loopType = $steps[$current]->{'stepLoop'};
		my $outputStep = $steps[$current]->{'output'};
		my $subref = \&$fname;
		if ($loopType eq 'sample') {
			for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
				my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
				my $outputLocation = $outputStep. "/" .$sampleName;
				SubmitToCluster::initSubmit(\%cfg, $outputLocation);
				# Tests for the first step in the list. Used for dependencies.
				my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, $sampleName, $rAoH_sampleLanes, $rAoH_seqDictionary, \%jobIdVarPrefix);
				if (defined($jobIdVar)) {
					$globalDep{$fname}->{$sampleName} = $jobIdVar;
				}
			}
		}elsif ($loopType eq 'design') {
			for my $design (keys %{$rHoAoA_designGroup}) {
					my $outputLocation = $outputStep. "/" .$design. '.0';
					SubmitToCluster::initSubmit(\%cfg, $outputLocation);
					my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, $design, $rHoAoA_designGroup, $rHoAoH_sampleInfo , \%jobIdVarPrefix);
					if (defined($jobIdVar)) {
						$globalDep{$fname}->{$design} = $jobIdVar;
					}
			}
		}else {
			my $outputLocation = $outputStep;
			SubmitToCluster::initSubmit(\%cfg, $outputLocation);
			# Tests for the first step in the list. Used for dependencies.
			my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, $designFilePath, $rHoAoA_designGroup, $rHoAoH_sampleInfo , \%jobIdVarPrefix);
			if (defined($jobIdVar)) {
				$globalDep{$fname}->{$fname} = $jobIdVar;
			}
		}
	}  
}

sub trimming {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;
	my $trimJobIdVarNameSample = undef;
	my $inputFile;
	my $outputFile;
	my $outputFastqPair1Name;
	
	for my $rH_laneInfo (@$rAoH_sampleLanes) {
			
		my $minQuality  = $rH_cfg->{'trim.minQuality'};
		my $minLength   = $rH_cfg->{'trim.minLength'};
		my $laneDirectory = 'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
		print "mkdir -p $laneDirectory\n";
		
		my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo, $laneDirectory);
		my $trimJobIdVarNameLane=undef;
		if( $rH_trimDetails->{'command'} ne "") {
			$trimJobIdVarNameLane = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $rH_trimDetails->{'command'}, 'reads/' .$sampleName, $workDirectory);
			$trimJobIdVarNameLane = '$' .$trimJobIdVarNameLane ;
			$trimJobIdVarNameSample .= $trimJobIdVarNameLane .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	}
	if (defined($trimJobIdVarNameSample )){
	 $trimJobIdVarNameSample = substr $trimJobIdVarNameSample, 0, -1 ;
	}
	return $trimJobIdVarNameSample; 
}


sub rawCounts {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;
  my $alignJobIdVarNameSample = undef;
	my $jobDependency = undef;
	
	if($depends > 0) {
					$jobDependency = $globalDep{'aligning'}{$sampleName};
	}
	print "mkdir -p raw_counts/$sampleName/output_jobs\n";
	## Generate aligment stats
	for my $rH_laneInfo (@$rAoH_sampleLanes) {		
# 			##get raw read count	
# 		if ((defined(LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir'))) && (LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir') ne '')){
# 			$inputFile = LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir') .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};
# 		}else{
# 			$inputFile = $sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};
# 		}
# 		$outputFile= 'raw_counts/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.raw.csv' ;
# 		my $command = Metrics::readStats($rH_cfg,$inputFile,$outputFile,'fastq');
# 		my $rawReadStatJobID = undef;
# 		if( $command ne "") {
# 			$rawReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "raw_counts", undef, 'RAWREADSTAT' .$rH_jobIdPrefixe ->{$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $command, 'raw_counts/' .$sampleName, $workDirectory);
# 			$rawReadStatJobID = '$'.$rawReadStatJobID;
# 		}	
		my $command = undef;
		my $inputFile = undef;
		my $outputFile = undef;
		my $filteredReadStatJobID = undef;
		my $alignedReadStatJobID = undef;
		## get trimmed read count
		$inputFile = 'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.trim.out';
		$outputFile= 'raw_counts/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.filtered.csv' ;
		$command = Metrics::readStats($rH_cfg,$inputFile,$outputFile, $sampleName,'trim');
		
		if( defined($command) ) {
		  if ($command ne ""){
			$filteredReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "rawCounts", "filtered" , 'FILTERREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} ,$jobDependency , $sampleName, $command, 'raw_counts/'  .$sampleName, $workDirectory);
			$filteredReadStatJobID = '$'.$filteredReadStatJobID;			
			}
		}
		# get aligned read counts
		$inputFile = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam';
		$outputFile= 'raw_counts/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.aligned.csv' ;
		$command = Metrics::readStats($rH_cfg,$inputFile,$outputFile, $sampleName,'bam');		
		if( defined($command)) {
			if( $command ne "") {
			$alignedReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "rawCounts", "aligned", 'ALIGNEDREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , $jobDependency , $sampleName, $command, 'raw_counts/' .$sampleName, $workDirectory);
			$alignedReadStatJobID = '$'.$alignedReadStatJobID;
			}	
		}
		## Merge read stats
		#my $rawFile = 'raw_counts/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.raw.csv' ;
		my $rawFile = "";
		my $filterFile = 'raw_counts/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.filtered.csv' ;
		my $alignFile = $outputFile ;
		my $sampleNameFull = $sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
		$outputFile= 'raw_counts/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.csv' ;
		$command = Metrics::mergeIndvidualReadStats($rH_cfg, $sampleNameFull, $rawFile, $filterFile, $alignFile, $outputFile);
		my $mergeReadStatJobID = undef;
		if( defined($command) ) {
		  if (defined ($alignedReadStatJobID) && defined ($filteredReadStatJobID) && $command ne "" ){
			$mergeReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeCounts", undef, 'MERGEREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} ,$alignedReadStatJobID.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$filteredReadStatJobID, $sampleName, $command, 'raw_counts/' .$sampleName, $workDirectory);
			$mergeReadStatJobID = '$'.$mergeReadStatJobID;
			$alignJobIdVarNameSample .= $mergeReadStatJobID .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
			}
		}
	}
	if(defined($alignJobIdVarNameSample)){
		$alignJobIdVarNameSample = substr $alignJobIdVarNameSample, 0, -1 ;
	}
	return $alignJobIdVarNameSample;
}

sub aligning{
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefixe = shift;
 	my $jobDependency = undef;
 	my $bwaJobId;
 	
	if($depends > 0) {
		$jobDependency = $globalDep{'trimming'}{$sampleName};
	}
  print "BWA_JOB_IDS=\"\"\n";
	print "mkdir -p $sampleName/output_jobs reads/$sampleName/output_jobs raw_counts/$sampleName/output_jobs alignment/$sampleName/output_jobs\n";
	for my $rH_laneInfo (@$rAoH_sampleLanes) {
		my $alignJobIdVarNameLane;
		my $pair1;
		my $pair2;
		my $single;
		my $command;		
    my $rgId = $rH_laneInfo->{'libraryBarcode'} . "_" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rgSampleName = $rH_laneInfo->{'name'};
    my $rgLibrary = $rH_laneInfo->{'libraryBarcode'};
    my $rgPlatformUnit = 'run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rgCenter = LoadConfig::getParam( $rH_cfg, 'aln', 'bwaInstitution' );

		#align lanes
    my $outputAlnDir = 'alignment/'.$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    print 'mkdir -p '.$outputAlnDir."\n";
    my $outputAlnPrefix = $outputAlnDir.'/'.$sampleName;
    my $bwaJobId = undef ; 

		if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
			$single =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.single.fastq.gz';
      $pair1  = undef;
      $pair2  = undef;

	  }elsif($rH_laneInfo->{'runType'} eq "PAIRED_END") {
	    $single =  undef;
			$pair1 =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.pair1.fastq.gz';
			$pair2 =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.pair2.fastq.gz';	  	  
	  }
		
    my $rA_commands = BWA::aln($rH_cfg, $sampleName, $pair1, $pair2, $single, $outputAlnPrefix, $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter);
    if(@{$rA_commands} == 3) {
      my $read1JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read1.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ1ALN', $jobDependency, $sampleName, $rA_commands->[0], 'alignment/'. $sampleName , $workDirectory);
      $read1JobId = '$'.$read1JobId;
      my $read2JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read2.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ2ALN', $jobDependency, $sampleName, $rA_commands->[1], 'alignment/'. $sampleName , $workDirectory );
      $read2JobId = '$'.$read2JobId;
      $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'sampe.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA', $read1JobId.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$read2JobId, $sampleName, $rA_commands->[2], 'alignment/'. $sampleName , $workDirectory );
      $bwaJobId = '$'.$bwaJobId;
      print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$bwaJobId."\n";
    }
    else {
      my $readJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READALN', $jobDependency, $sampleName, $rA_commands->[0], 'alignment/'. $sampleName );
      $readJobId = '$'.$readJobId;
      $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'samse.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA',  $readJobId, $sampleName, $rA_commands->[1], 'alignment/'. $sampleName , $workDirectory);
      $bwaJobId = '$'.$bwaJobId;
      print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$bwaJobId."\n";
    } 
	}	
	# Merge /sort reads
	$jobDependency = '$BWA_JOB_IDS';
	my $latestBam;
  my @inputBams;
  my $outputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam';
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $directory = 'alignment/'.$sampleName."/run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'}."/";
    my $sortedLaneBamFile = $directory.$rH_laneInfo->{'name'}.".sorted.bam";
    my $runName = $sampleName."_run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'};
    push(@inputBams, $sortedLaneBamFile);
  }
  
  my $command = Picard::mergeFiles($rH_cfg, $sampleName, \@inputBams, $outputBAM);
  my $mergeJobId = undef;
  if( defined($command) && $command ne "" ) {
    $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeLanes", undef, 'MERGELANES'.$rH_jobIdPrefixe ->{$sampleName}, $jobDependency, $sampleName, $command, $workDirectory.'/'.$sampleName );
    $mergeJobId  = '$'.$mergeJobId;
  }
  return $mergeJobId;
}
    
sub wiggleRNASEQ {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'qcTagDirectories'}{$sampleName};
	}

	my $inputBAM       =  'alignment/' .$sampleName . '/' . $sampleName . '.sorted.bam' ; 
	my $outputBAM      =  $inputBAM;
	my $outputBedGraph =  'tracks/' . $sampleName . '/' . $sampleName . '.bedGraph';
	my $outputWiggle   =  'tracks/' . $sampleName . '/' . $sampleName . '.bw';
	my $prefixJobName  =  undef;
	
	print "mkdir -p tracks/$sampleName/output_jobs\n";

	my $wiggleJobId = undef ;
	my $command = Wiggle::graph($rH_cfg, $sampleName, $inputBAM, $outputBedGraph, $outputWiggle);
	if( $command ne "" ) {
		$wiggleJobId  = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", $prefixJobName, 'WIGGLE' .$rH_jobIdPrefixe ->{$sampleName} ,$jobDependency, $sampleName, $command, 'tracks/' .$sampleName, $workDirectory);
		$wiggleJobId  = '$'.$wiggleJobId ;
	}
	return $wiggleJobId;    
}


sub qcTagDirectories {
	my $depends     = shift;
	my $rH_cfg      = shift;
	my $sampleName  = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	my $qcTagsJobID   = undef;

	if($depends > 0) {
		$jobDependency = $globalDep{'aligning'}{$sampleName};
	}

	my $inputBAM = 'alignment/' . $sampleName . "/" . $sampleName . ".sorted.bam" ; 
	
	# Create command
	my $rA_commands = HOMER::makeTagDirectory($rH_cfg, $sampleName, $inputBAM , 'tags');
	if(@{$rA_commands} == 3) {
		my $qcTagsJobID1  = SubmitToCluster::printSubmitCmd($rH_cfg, "qcTags", 1, 'QCTAGSTMP' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency , $sampleName, $rA_commands->[0], 'tags/' .$sampleName, $workDirectory);
		$qcTagsJobID1  = '$'.$qcTagsJobID1 ;
 		my $qcTagsJobID2  = SubmitToCluster::printSubmitCmd($rH_cfg, "qcTags", 2, 'QCTAGSFILLTMP' .$rH_jobIdPrefixe ->{$sampleName} , $qcTagsJobID1 , $sampleName, $rA_commands->[1], 'tags/' .$sampleName, $workDirectory);
		$qcTagsJobID2  = '$'.$qcTagsJobID2 ;
 		my $qcTagsJobID3  = SubmitToCluster::printSubmitCmd($rH_cfg, "qcTags", 3, 'QCTAGSCREATETAGS' .$rH_jobIdPrefixe ->{$sampleName} , $qcTagsJobID1 , $sampleName, $rA_commands->[2], 'tags/' .$sampleName, $workDirectory);
		$qcTagsJobID3  = '$'.$qcTagsJobID3 ;
		$qcTagsJobID   = $qcTagsJobID3;
	}elsif(@{$rA_commands} == 1){
		my $qcTagsJobID   = SubmitToCluster::printSubmitCmd($rH_cfg, "qcTags", undef, 'QCTAGS' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency , $sampleName, $rA_commands->[0], 'tags/' .$sampleName, $workDirectory);
		$qcTagsJobID   = '$'.$qcTagsJobID ;
	}
	print 'QCTAGS_JOB_IDS=${QCTAGS_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$qcTagsJobID."\n";
	return $qcTagsJobID;    
}
sub wiggle {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'qcTagDirectories'}{$sampleName};
	}

	my $tagDirectory   =  'tags/'. $sampleName;
	my $outputWiggle   =  'tracks/' . $sampleName . '/' . $sampleName . '.ucsc.bedGraph.gz';
	my $prefixJobName  =  undef;
	
	print "mkdir -p tracks/$sampleName/output_jobs\n";

	my $wiggleJobId = undef ;
	my $command = HOMER::makeUCSCFile($rH_cfg, $sampleName, $tagDirectory , $outputWiggle);
	if(defined($command) && $command ne "") {
		$wiggleJobId  = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", $prefixJobName, 'WIGGLE' .$rH_jobIdPrefixe ->{$sampleName} ,$jobDependency, $sampleName, $command, 'tracks/' .$sampleName, $workDirectory);
		$wiggleJobId  = '$'.$wiggleJobId ;
	}
	return $wiggleJobId;    
}


sub qcPlots {
	my $depends = shift;
	my $rH_cfg = shift;
	my $designFilePath = shift;
	
	my $jobDependency = undef;
	if($depends > 0) {
			$jobDependency = '$QCTAGS_JOB_IDS';
	}
	print "mkdir -p graphs/output_jobs\n";
	my $qcPlotsJobId = undef ;
	
	my $command = HOMER::qcPlotsR($rH_cfg, $designFilePath, $workDirectory);
	if(defined($command) && $command ne "") {
		$qcPlotsJobId  = SubmitToCluster::printSubmitCmd($rH_cfg, "qcGraphs", undef, 'QCPLOTS' , $jobDependency, undef, $command, 'graphs/', $workDirectory);
		$qcPlotsJobId  = '$'.$qcPlotsJobId ;
	}
	return $qcPlotsJobId;    

}

sub peakCall {
	my $depends = shift;
	my $rH_cfg = shift;
	my $design = shift;
	my $rHoAoA_designGroup  = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rH_jobIdPrefixe = shift;
	
	my $jobDependency = undef;
  my $PeakCallJobIdVarNameSample = undef;
  
	## Design
	# Process Treatments ( and controls when applicable)
	my $numberTreatments 	=  @{$rHoAoA_designGroup->{$design}->[0]};	
	my $numberControls 		=  @{$rHoAoA_designGroup->{$design}->[1]};	
	my $control;
	my $treatment;
	my $controlBam;
	my $treatmentBam;
	
	if ($numberTreatments >= 1) {
		# At least one treatment
		for (my $j = 0;   $j < $numberTreatments; $j++) {
			if($numberControls == 1){
					$control       = $rHoAoA_designGroup->{$design}->[1]->[0];
					$controlBam    = 'alignment/' .$control . '/' . $control . '.sorted.bam' ;
			}elsif( $numberControls == $numberTreatments){
					$control       = $rHoAoA_designGroup->{$design}->[1]->[$j];
					$controlBam    = 'alignment/' .$control . '/' .$control . '.sorted.bam' ; 
			}else{
					$control       = undef ;
					$controlBam    = undef;
			}
			$treatment    =  $rHoAoA_designGroup->{$design}->[0]->[$j];
			$treatmentBam =  'alignment/' .$treatment . '/' . $treatment . '.sorted.bam' ;
			
			my $outputPath = 'peak_call/' .$design. '.' .$j;
			if(defined($control)){
				print "# design ".$design. ", treatment= ".$treatment.", control=". $control. "\n" ;
			}else{
				print "# design ".$design. ", treatment= ".$treatment."". "\n" ;
			}
			# Run type for treatment
			my $paired = undef ; 
			my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$treatment};
			for my $rH_laneInfo (@$rAoH_sampleLanes) {
				if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
					$paired = 0;
				}elsif( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
			    $paired = 1;
				}else{
					die "Unknown runType: " . $rH_laneInfo->{'runType'} . "\n";
				}
				last if (defined($paired));
			}
			if(! defined($rHoAoH_sampleInfo->{$treatment}) ){
				die "ERROR: Sample ".$treatment ." in design file is not present in NANUQ sample file". "\n";
			}elsif(defined($control) && !defined($rHoAoH_sampleInfo->{$control})){
				die "ERROR: Sample ".$control ." in design file is not present in NANUQ sample file ". "\n";
			}
			# MACS command
			my $command = MACS2::generatePeaks( $rH_cfg, $design, $treatmentBam, $controlBam , $rHoAoA_designGroup->{$design}->[2]->[0], $outputPath, $paired );
			if( defined($command) && $command ne "") {
			  if( $depends > 0 && defined( $control ) && defined($globalDep{'aligning'}{ $control }) && defined($globalDep{'aligning'}{ $treatment })){
					$jobDependency = $globalDep{'aligning'}{ $control }.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$globalDep{'aligning'}{ $treatment };
				}elsif( $depends > 0 && defined($globalDep{'aligning'}{ $treatment })) {
				  $jobDependency = $globalDep{'aligning'}{ $treatment };
				}
				print "mkdir -p " . $outputPath . '/output_jobs/' . "\n";
				my $peakCallJobId= SubmitToCluster::printSubmitCmd($rH_cfg, 'peakCall', $j, 'MACS2PEAKCALL'. $rH_jobIdPrefixe ->{$design}. $j , $jobDependency, $design, $command, $outputPath , $workDirectory);
				$peakCallJobId= '$' .$peakCallJobId;
				$PeakCallJobIdVarNameSample .= $peakCallJobId.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
			}
		}
	} else{
		; # do nothing
	}
	if(defined($PeakCallJobIdVarNameSample )){
	  $PeakCallJobIdVarNameSample = substr $PeakCallJobIdVarNameSample, 0, -1 ;
	}
	return $PeakCallJobIdVarNameSample;
}

sub annotation{
	my $depends = shift;
	my $rH_cfg = shift;
	my $design = shift;
	my $rHoAoA_designGroup  = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'peakCall'}{$design};
	}	
	##Iterate over design
	my $numberTreatments 	=  @{$rHoAoA_designGroup->{$design}->[0]};	
	my $annotationJobIdGroups = undef;
	if ($numberTreatments >= 1) {
		# At least one treatment
		for (my $j = 0;   $j < $numberTreatments; $j++) {
		  my $peakCallPath = 'peak_call/' .$design. '.' .$j ;
			my $outputPath = 'annotation/' .$design. '.' .$j ;
			my $peakType 	=  $rHoAoA_designGroup->{$design}->[2]->[0];
			if ($peakType eq "N"){
				my $bedName = $peakCallPath . '/' .$design. '_peaks.bed';
				my $command = HOMER::annotatePeaks($rH_cfg, $design, $bedName, $outputPath);
				if( defined($command) && $command ne "") {
					print "mkdir -p " . $outputPath . '/output_jobs/' . "\n";
					my $annotationJobId= SubmitToCluster::printSubmitCmd($rH_cfg, 'annotation', undef, 'HOMERANNOTATION' .$rH_jobIdPrefixe->{$design} , $jobDependency, $design, $command, $outputPath , $workDirectory);
					$annotationJobId= '$' .$annotationJobId;
					$annotationJobIdGroups.= $annotationJobId.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
				}
			}else{
				; # do nothing
			}
		}
	}
	if(defined($annotationJobIdGroups)){
	  $annotationJobIdGroups = substr $annotationJobIdGroups, 0, -1 ;
	}
	return $annotationJobIdGroups;
} 
# sub printSubmitCmd $rH_cfg	$stepName       $jobNameSuffix  $jobIdPrefix    $dependancyName $sampleName  $command  $outputDir $workDirectory 
sub motif{
	my $depends = shift;
	my $rH_cfg = shift;
	my $design = shift;
	my $rHoAoA_designGroup  = shift;
	my $rHoAoH_sampleInfo  = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'peakCall'}{$design};
	}	
		##Iterate over design
	my $numberTreatments 	=  @{$rHoAoA_designGroup->{$design}->[0]};	
	my $motifJobIdGroups=undef;
	if ($numberTreatments >= 1) {
		# At least one treatment
		for (my $j = 0;   $j < $numberTreatments; $j++) {	
			my $outputPath = 'motif/' .$design. '.' .$j ;
			my $peakCallPath = 'peak_call/' .$design. '.' .$j ;
			my $peakType 	=  $rHoAoA_designGroup->{$design}->[2]->[0];
			if ($peakType eq "N"){
				my $bedName = $peakCallPath . '/' .$design. '_peaks.bed';
				my $command = HOMER::generateMotif($rH_cfg, $design, $bedName, $outputPath);
				if( defined($command) && length($command) > 0) {
				  print "mkdir -p " . $outputPath . '/output_jobs/' . "\n";
					my $motifJobId= SubmitToCluster::printSubmitCmd($rH_cfg, 'motif', undef, 'HOMERMOTIF' .$rH_jobIdPrefixe->{$design} , $jobDependency, $design, $command, $outputPath , $workDirectory);
					$motifJobId= '$' .$motifJobId;
					$motifJobIdGroups.= $motifJobId.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
				}
			}else{
				; # do nothing
			}	
		}
	}			
	if(defined($motifJobIdGroups)){
	  $motifJobIdGroups = substr $motifJobIdGroups, 0, -1 ;
	}
	return $motifJobIdGroups;
}

1;