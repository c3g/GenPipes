#!/usr/bin/perl

=head1 NAME

I<rnaSeq>

=head1 SYNOPSIS

rnaSeq.pl

=head1 DESCRIPTION

B<rnaSeq> Is the main RNAseq pipeline.

=head1 AUTHOR

B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

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

}


# Dependencies
#--------------------
use Getopt::Std;

use LoadConfig;
use Picard;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SubmitToCluster;
use TophatBowtie;
use Trimmomatic;
use Metrics;
use Cufflinks;
use Wiggle;
use HtseqCount;
use DiffExpression;

#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'trimming' , 'stepLoop' => 'sample' , 'output' => 'reads'});
push(@steps, {'name' => 'aligning' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
push(@steps, {'name' => 'merging' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
#push(@steps, {'name' => 'mutation' , 'stepLoop' => 'sample' , 'output' => 'mpileup'});
push(@steps, {'name' => 'wiggle' , 'stepLoop' => 'sample' , 'output' => 'tracks'});
push(@steps, {'name' => 'rawCounts' , 'stepLoop' => 'sample' , 'output' => 'raw_count'});
push(@steps, {'name' => 'fpkm' , 'stepLoop' => 'sample' , 'output' => 'fpkm'});
push(@steps, {'name' => 'cuffdiff' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'metrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
push(@steps, {'name' => 'dge' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'goseq' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'delivrable' , 'stepLoop' => 'group' , 'output' => 'Delivrable'});


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
	##get design groups
	my $rHoAoA_designGroup = Cufflinks::getDesign(\%cfg,$designFilePath);
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
	#generate design jobIdprefix
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
		}
		else {
			my $outputLocation = $outputStep;
			SubmitToCluster::initSubmit(\%cfg, $outputLocation);
			# Tests for the first step in the list. Used for dependencies.
			my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, $rHoAoH_sampleInfo, $rHoAoA_designGroup, $rAoH_seqDictionary, \%jobIdVarPrefix);
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
	my $libraryType = LoadConfig::getParam($rH_cfg, 'default', 'libraryType');
	for my $rH_laneInfo (@$rAoH_sampleLanes) {
		print "mkdir -p metrics/$sampleName/output_jobs reads/$sampleName/output_jobs\n";
		##get raw read count
# 		my $inputFile = LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir') .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};
# 		my $outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.raw.csv' ;
# 		my $command = Metrics::readStats($rH_cfg,$inputFile,$outputFile,'fastq',$libraryType);
# 		my $rawReadStatJobID = undef;
# 		if(defined($command) && length($command) > 0) {
# 			$rawReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'raw', 'RAWREADSTAT' .$rH_jobIdPrefixe ->{$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $command, 'metrics/' .$sampleName, $workDirectory);
# 			$rawReadStatJobID = '$'.$rawReadStatJobID;
# 		}
# 		
		## trimming - TO DO should be modified to the new rawread location (cf. David modif) and portability
		my $minQuality  = $rH_cfg->{'trim.minQuality'};
		my $minLength   = $rH_cfg->{'trim.minLength'};
		my $laneDirectory = 'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
		print "mkdir -p $laneDirectory\n";
		my $outputFastqPair1Name;
		if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
			$outputFastqPair1Name = $laneDirectory . $sampleName.'.t'.$minQuality.'l'.$minLength.'.single.fastq.gz';
		}
		elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
			$outputFastqPair1Name = $laneDirectory . $sampleName.'.t'.$minQuality.'l'.$minLength.'.pair1.fastq.gz';
		}
		else {
			die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
		}
		my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo, $laneDirectory);
		my $trimJobIdVarNameLane=undef;
		if(length($rH_trimDetails->{'command'}) > 0) {
			$trimJobIdVarNameLane = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $rH_trimDetails->{'command'}, 'reads/' .$sampleName, $workDirectory);
			$trimJobIdVarNameLane = '$' .$trimJobIdVarNameLane ;
		}
		my $trinityOut = $laneDirectory .'/' . $sampleName . '.trim.out';
		##get trimmed read count
		my $outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.triming.tsv' ;
		my $command = Metrics::readStats($rH_cfg,$trinityOut,$outputFile,$sampleName,'trim');
		my $filteredReadStatJobID ;
		if(defined($command) && length($command) > 0) {
			$filteredReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'filtered', 'FILTERREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} ,$trimJobIdVarNameLane, $sampleName, $command,'metrics/'  .$sampleName, $workDirectory);
			$filteredReadStatJobID = '$'.$filteredReadStatJobID;
			$trimJobIdVarNameSample .= $filteredReadStatJobID .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	}
	if(defined($trimJobIdVarNameSample) && length($trimJobIdVarNameSample) > 0) {
		$trimJobIdVarNameSample = substr $trimJobIdVarNameSample, 0, -1 ;
	}
	return $trimJobIdVarNameSample;	
}

sub aligning {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $alignJobIdVarNameSample = undef;
	my $jobDependency = undef;
	if($depends > 0) {
	$jobDependency = $globalDep{'trimming'}{$sampleName};
	}
	
	print "mkdir -p reads/$sampleName/output_jobs metrics/$sampleName/output_jobs alignment/$sampleName/output_jobs\n";
	for my $rH_laneInfo (@$rAoH_sampleLanes) {
		my $alignJobIdVarNameLane;
		my $pair1;
		my $pair2;
		my $single;
		my $command;
		#align lanes
		my $workDir = 'reads' ;
		my $outputDirPath = 'alignment/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
		print "mkdir -p $outputDirPath \n" ;
		if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
			$single =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.single.fastq.gz';
			$command = TophatBowtie::align($rH_cfg, $sampleName, $rH_laneInfo, $single, ' ' );
		}
		elsif($rH_laneInfo->{'runType'} eq "PAIRED_END") {
			$pair1 =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.pair1.fastq.gz';
			$pair2 =  'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName .'.t' .LoadConfig::getParam($rH_cfg,'trim','minQuality') .'l' .LoadConfig::getParam($rH_cfg,'trim','minLength') .'.pair2.fastq.gz';
			$command = TophatBowtie::align($rH_cfg, $sampleName, $rH_laneInfo, $pair1, $pair2);
		}
		if(defined $command && length($command) > 0){
			$alignJobIdVarNameLane = SubmitToCluster::printSubmitCmd($rH_cfg, "align", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'ALIGN' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} .'ALIGN', $jobDependency, $sampleName, $command, 'alignment/' .$sampleName, $workDirectory);
			$alignJobIdVarNameSample .= '$'. $alignJobIdVarNameLane .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'); 
		} 
		##generate aligment stats
# 		my $inputFile = 'alignment/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . 'accepted_hits.bam';
# 		my $outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.aligned.csv' ;
# 		$command = Metrics::readStats($rH_cfg,$inputFile,$outputFile,'bam');
# 		my $alignedReadStatJobID = undef;
# 		if(defined($command) && length($command) > 0) {
# 			$alignedReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'aligned', 'ALIGNEDREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} ,$alignJobIdVarNameLane, $sampleName, $command, 'metrics/' .$sampleName, $workDirectory);
# 			$alignedReadStatJobID = '$'.$alignedReadStatJobID;
# 		}
# 		##merge read stats
# 		my $rawFile = 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.raw.csv' ;
# 		my $filterFile = 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.filtered.csv' ;
# 		my $alignFile = $outputFile ;
# 		my $sampleNameFull = $sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
# 		$outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.csv' ;
# 		$command = Metrics::mergeIndvidualReadStats($rH_cfg, $sampleNameFull, $rawFile, $filterFile, $alignFile, $outputFile);
# 		my $mergeReadStatJobID = undef;
# 		if(defined($command) && length($command) > 0) {
# 			$mergeReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'merged', 'MERGEREADSTAT' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} ,$alignedReadStatJobID, $sampleName, $command, 'metrics/' .$sampleName, $workDirectory);
# 			$mergeReadStatJobID = '$'.$mergeReadStatJobID;
# 			$alignJobIdVarNameSample .= $mergeReadStatJobID .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
# 		}
		
	}
	$alignJobIdVarNameSample = substr $alignJobIdVarNameSample, 0, -1 ;
	return $alignJobIdVarNameSample;
}

sub merging {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'aligning'}{$sampleName};
	}
	##Merging
	my $inputBAM ; 
	my $outputBAM = 'alignment/' .$sampleName .'/' .$sampleName .'.merged.bam' ;
	my $workDir = 'alignment' ;
	my @alignFiles;
	print "mkdir -p alignment/$sampleName/output_jobs\n";
	for my $rH_laneInfo (@$rAoH_sampleLanes) {
		my $laneDirectory = "alignment/" . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
		$inputBAM = $laneDirectory . 'accepted_hits.bam';
		push(@alignFiles, $inputBAM) ;
	}
	my $command = Picard::mergeFiles($rH_cfg, $sampleName, \@alignFiles, $outputBAM);
	my $mergeJobId = undef;
	if(defined($command) && length($command) > 0) {
		$mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeFiles", undef, 'MERGELANES' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $command, $workDir .'/' .$sampleName, $workDirectory);
		$mergeJobId = '$'.$mergeJobId;
	}
	## reorder
	$inputBAM = $outputBAM;
	$outputBAM = 'alignment/' .$sampleName .'/' .$sampleName .'.merged.karyotypic.bam';
	$command = Picard::reorderSam($rH_cfg, $sampleName, $inputBAM, $outputBAM);
	my $reorderJobId = undef;
	if(defined($command) && length($command) > 0) {
		$reorderJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "reorderSam", undef, 'REORDER' .$rH_jobIdPrefixe ->{$sampleName} .'REORDER', $mergeJobId, $sampleName, $command, $workDir .'/' .$sampleName, $workDirectory);
		$reorderJobId = '$'.$reorderJobId;
	}
	## mark duplicates
	$inputBAM = $outputBAM;
	$outputBAM = 'alignment/' .$sampleName .'/' .$sampleName .'.merged.mdup.bam';
	my $duplicatesMetricsFile = 'alignment/' .$sampleName .'/' .$sampleName .'.merged.mdup.metrics';
	$command = Picard::markDup($rH_cfg, $sampleName, $inputBAM, $outputBAM,$duplicatesMetricsFile );
	my $markDupJobId = undef;
	if(defined($command) && length($command) > 0) {
		$markDupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'MARKDUP' .$rH_jobIdPrefixe ->{$sampleName} , $reorderJobId, $sampleName, $command, $workDir .'/' .$sampleName, $workDirectory);
		$markDupJobId = '$'.$markDupJobId
	}
	return $markDupJobId;
}

# sub mutation{
# 	my $depends = shift;
# 	my $rH_cfg = shift;
# 	my $sampleName = shift;
# 	my $rAoH_sampleLanes  = shift;
# 	my $rAoH_seqDictionary = shift;
# ALIGNEDREADSTAT_JOB_ID
# 	my $jobDependency = undef;
# 	if($depends > 0) {
# 		$jobDependency = $globalDep{'merging'}{$sampleName};
# 	}
# 
# 	my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.bam" ; 
# 
# 
# }

sub wiggle {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}

	my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.bam" ; 
	#testing for strand-specificity

	my $strandSPecificityInfo = LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
	my @strandJobId ;
	my @outputBAM;
	my @outputBedGraph;
	my @outputWiggle;
	my @prefixJobName;
	print "mkdir -p alignment/$sampleName/output_jobs tracks/$sampleName/output_jobs\n";
	if($strandSPecificityInfo ne "fr-unstranded") {
	## strand specific 
		@outputBAM = {'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.forward.bam' ,  'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.reverse.bam'};
		my $rA_command = Wiggle::strandBam($rH_cfg, $sampleName, $inputBAM, \@outputBAM);
		if(defined($rA_command) && @{$rA_command} > 1) {
			my $strandJobIdF = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", 'FORWARD', 'FSTRANDSPEC' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $rA_command->[0], 'alignment/' .$sampleName, $workDirectory);
			push(@strandJobId, '$'.$strandJobIdF );
			my $strandJobIdR = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", 'REVERSE', 'RSTRANDSPEC' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $rA_command->[1], 'alignment/' .$sampleName, $workDirectory);
			push(@strandJobId, '$'.$strandJobIdR );
		}
		@outputBedGraph = {'tracks/' . $sampleName . '/' . $sampleName . '.forward.bedGraph' ,  'tracks/' . $sampleName . '/' . $sampleName . '.reverse.bedGraph'};
		@outputWiggle = {'tracks/' . $sampleName . '/' . $sampleName . '.forward.bw' ,  'tracks/' . $sampleName . '/' . $sampleName . '.reverse.bw'};
		@prefixJobName = { 'FORWARD', 'REVERSE'};
	}
	else {
		push(@outputBAM,$inputBAM);
		push(@strandJobId, $jobDependency);
		push(@outputBedGraph,'tracks/' . $sampleName . '/' . $sampleName . '.bedGraph');
		push(@outputWiggle,'tracks/' . $sampleName . '/' . $sampleName . '.bw' );
		push(@prefixJobName , undef ) ;
	}
	my $wiggleJobId ;
	for(my $i = 0; $i <@outputBAM; $i++) {
		my $command = Wiggle::graph($rH_cfg, $sampleName, $inputBAM, $outputBedGraph[$i], $outputWiggle[$i]);
		if(defined($command) && length($command) > 0) {
			my $tmpJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", $prefixJobName[$i], 'WIGGLE' .$rH_jobIdPrefixe ->{$sampleName} , $strandJobId[$i], $sampleName, $command, 'tracks/' .$sampleName, $workDirectory);
			$wiggleJobId .= '$'.$tmpJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	} 
	$wiggleJobId = substr $wiggleJobId, 0, -1 ;
	return $wiggleJobId;	
}


sub rawCounts {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}
	print "mkdir -p alignment/$sampleName/output_jobs raw_counts/$sampleName/output_jobs\n";
	my $inputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam' ;
        my $sortedBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.queryNameSorted.bam' ;
	my $inputGtf = LoadConfig::getParam($rH_cfg, 'htseq', 'referenceGtf');
	my $outputCount = 'raw_counts/' . $sampleName . '.readcounts.csv';
	my $sortOrder = 'queryname';
	my $strandInfo;
	my $strandSPecificityInfo = LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
	if($strandSPecificityInfo ne "fr-unstranded") {
		 $strandInfo= 'yes';
	}
	else {
		$strandInfo= 'no';
	}
	## query sort the bam
	my $sortJobId;
	my $command = Picard::sortSam($rH_cfg, $sampleName, $inputBAM, $sortedBAM, $sortOrder);
	if(defined($command) && length($command) > 0) {
		$sortJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "sortSam", undef, 'QNSORT' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $command, 'alignment/' .$sampleName, $workDirectory);
		$sortJobId='$'.$sortJobId
	}
	## count reads
        my $countJobId;
	$command = HtseqCount::readCountPortable($rH_cfg, $sortedBAM, $inputGtf, $outputCount, $strandInfo); 
	if(defined($command) && length($command) > 0) {
		$countJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "htseq", undef, 'RAWCOUNT' .$rH_jobIdPrefixe ->{$sampleName} , $sortJobId, $sampleName, $command, 'raw_counts/' .$sampleName, $workDirectory);
		$countJobId='$'.$countJobId
	}
	return $countJobId;
}


sub fpkm {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}
	print "mkdir -p fpkm/known/$sampleName fpkm/denovo/$sampleName fpkm/$sampleName/output_jobs\n";
	my $inputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam' ;
	my $outputKnown = 'fpkm/known/' . $sampleName;
	my $outputDeNovo = 'fpkm/denovo/' . $sampleName;
	my $gtfOption = '-G ' .LoadConfig::getParam($rH_cfg, 'fpkm','referenceGtf');
	
	
	## known FPKM
	my $fpkmJobId = undef;
	my $command = Cufflinks::fpkm($rH_cfg, $inputBAM, $outputKnown, $gtfOption);
	if(defined($command) && length($command) > 0) {
		my $fpkmKnownJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "fpkm", "KNOWN", 'FPKMK' .$rH_jobIdPrefixe ->{$sampleName} .'FPKM', $jobDependency, $sampleName, $command, 'fpkm/' .$sampleName, $workDirectory);
		$fpkmJobId = '$' .$fpkmKnownJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	## denovo FPKM
	$command = Cufflinks::fpkm($rH_cfg, $inputBAM, $outputDeNovo, undef);
	if(defined($command) && length($command) > 0) {
		my $fpkmDeNovoJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "fpkm", "DENOVO", 'FPKMD' .$rH_jobIdPrefixe ->{$sampleName} , $jobDependency, $sampleName, $command, 'fpkm/' .$sampleName, $workDirectory);
		$fpkmJobId .= '$' .$fpkmDeNovoJobId  .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	if (defined($fpkmJobId)) {
		$fpkmJobId = substr $fpkmJobId, 0, -1 ;
	}
	return $fpkmJobId;
}


sub cuffdiff {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;
	
	my $jobDependency = undef;
	if($depends > 0 and values(%{$globalDep{'fpkm'}}) > 0) {
		$jobDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'fpkm'}}));
	}
	print "mkdir -p cuffdiff/known cuffdiff/denovo cuffdiff/output_jobs\n";
	##iterate over design
	my $cuffddiffJobId;
	for my $design (keys %{$rHoAoA_designGroup}) {
		mkdir  $workDirectory ;
		mkdir  $workDirectory .'/cuffdiff';
		mkdir  $workDirectory .'/cuffdiff/denovo/' ;
		mkdir  $workDirectory .'/cuffdiff/denovo/' .$design ;
		## create the list of deNovo gtf to merge
		print "mkdir -p cuffdiff/$design/output_jobs\n";
		my $mergeListFile = $workDirectory .'/cuffdiff/denovo/' .$design .'/gtfMerge.list';
		open(MERGEF, ">$mergeListFile") or  die ("Unable to open $mergeListFile for wrtting") ;
		my $numberGroups = @{$rHoAoA_designGroup->{$design}} ;
		##iterate over group
		my @groupInuptFiles;
		for (my $i = 0;   $i < $numberGroups; $i++) {
			##iterate over samples in the design
			my $numberSample =  @{$rHoAoA_designGroup->{$design}->[$i]};
			my $gtfFile = '';
			my $bamfile = ' ';
			for (my $j = 0;   $j < $numberSample; $j++) {
				$gtfFile = 'fpkm/denovo/' .$rHoAoA_designGroup->{$design}->[$i]->[$j] .'/transcripts.gtf' ;
				print MERGEF $gtfFile;
				print MERGEF "\n";
				$bamfile .= 'alignment/' .$rHoAoA_designGroup->{$design}->[$i]->[$j] . '/' .$rHoAoA_designGroup->{$design}->[$i]->[$j] . '.merged.mdup.bam' .',' ;
			}
			$bamfile = substr $bamfile, 0, -1 ;
			push(@groupInuptFiles,$bamfile);
		}
		close($mergeListFile);

		my $outputPathKnown = 'cuffdiff/known/' .$design;
		my $outputPathDeNovo = 'cuffdiff/denovo/' .$design;
		
		my $command = Cufflinks::cuffmerge($rH_cfg, $mergeListFile, $outputPathDeNovo);
		my $cuffmergeJobId ;
		if(defined($command) && length($command) > 0) {
			$cuffmergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "cuffmerge", "MERGE", 'GTFMERGE' .$rH_jobIdPrefixe ->{$design} , $jobDependency, $design, $command, 'cuffdiff/' .$design, $workDirectory);
			$cuffmergeJobId = '$' .$cuffmergeJobId 
		}
		
		my $gtfDnMerged = 'cuffdiff/denovo/' .$design .'/merged.gtf';
		my $gtfDnFormatMerged = 'cuffdiff/denovo/' .$design .'/formated.merged.gtf';
		$command = Cufflinks::mergeGtfFormat($rH_cfg, $gtfDnMerged, $gtfDnFormatMerged);
		my $formatJobId;
		if(defined($command) && length($command) > 0) {
			$formatJobId= SubmitToCluster::printSubmitCmd($rH_cfg, "default", "FORMAT", 'GTFFORMAT' .$rH_jobIdPrefixe ->{$design} , $cuffmergeJobId, $design, $command, 'cuffdiff/' .$design, $workDirectory);
			$formatJobId= '$' .$formatJobId
		}

		##cuffdiff known
		$command = Cufflinks::cuffdiff($rH_cfg,\@groupInuptFiles,$outputPathKnown,LoadConfig::getParam($rH_cfg, 'cuffdiff','referenceGtf'));
		if(defined($command) && length($command) > 0) {
			my $cuffdiffKnownJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "cuffdiff", "KNOWN",  'CUFFDIFFK' .$rH_jobIdPrefixe ->{$design} , $jobDependency, $design, $command, 'cuffdiff/' .$design, $workDirectory);
			$cuffddiffJobId .= '$' .$cuffdiffKnownJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
		
		##cuffdiff de novo
		$command = Cufflinks::cuffdiff($rH_cfg,\@groupInuptFiles,$outputPathDeNovo,$gtfDnFormatMerged);
		if(defined($command) && length($command) > 0) {
			my $cuffdiffKnownJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "cuffdiff", "DENOVO", 'CUFFDIFFD' .$rH_jobIdPrefixe ->{$design} , $formatJobId, $design, $command, 'cuffdiff/' .$design, $workDirectory);
			$cuffddiffJobId .= '$' .$cuffdiffKnownJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	}
	$cuffddiffJobId = substr $cuffddiffJobId, 0, -1 ;
	my $command = Cufflinks::mergeCuffdiffRes($rH_cfg,$designFilePath,'cuffdiff','fpkm');
	my $mergeCuffdiffResJobID;
	if(defined($command) && length($command) > 0) {
		$mergeCuffdiffResJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "default", "MERGE_RES", 'CUFF_MERGE_RES', $cuffddiffJobId, undef, $command, 'cuffdiff/', $workDirectory);
		$mergeCuffdiffResJobID = '$' .$mergeCuffdiffResJobID;
	}
	$command = Cufflinks::filterResults($rH_cfg,'cuffdiff/known/') ;
	my $filterCuffdiffResJobID;
	if(defined($command) && length($command) > 0) {
		my $filterKCuffdiffResJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "default", "FILTERK", 'KFILT_CUFFDIFF', $mergeCuffdiffResJobID, undef, $command, 'cuffdiff/', $workDirectory);
		$filterCuffdiffResJobID .= '$' .$filterKCuffdiffResJobID .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	$command = Cufflinks::filterResults($rH_cfg,'cuffdiff/denovo/') ;
	if(defined($command) && length($command) > 0) {
		my $filterDCuffdiffResJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "default", "FILTERD", 'DFILT_CUFFDIFF', $mergeCuffdiffResJobID, undef, $command, 'cuffdiff/', $workDirectory);
		$filterCuffdiffResJobID .= '$' .$filterDCuffdiffResJobID .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	$filterCuffdiffResJobID= substr $filterCuffdiffResJobID, 0, -1 ;
	return $filterCuffdiffResJobID;
}


sub metrics {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $mergingDependency = undef;
	if($depends > 0) {
		$mergingDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'merging'}}));
	}
	## RNAseQC metrics
	mkdir  $workDirectory .'/alignment' ;
	open(RNASAMPLE, ">alignment/rnaseqc.samples.txt") or  die ("Unable to open alignment/rnaseqc.samples.txt for wrtting") ;
	print RNASAMPLE "Sample\tBamFile\tNote\n";
	my $projectName = LoadConfig::getParam($rH_cfg, 'metricsRNA', 'projectName');
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
		print RNASAMPLE "$sampleName\talignment/$sampleName/$sampleName.merged.mdup.bam\t$projectName\n";
	}
	print "mkdir -p metrics/output_jobs\n";
	my $sampleList = 'alignment/rnaseqc.samples.txt';
	my $outputFolder = 'metrics/';
	my $command = Metrics::rnaQc($rH_cfg, $sampleList, $outputFolder);
	my $rnaqcJobId = undef;
	my $metricsJobId = undef;
	if(defined($command) && length($command) > 0) {
		$rnaqcJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "rnaQc", undef, 'METRICSRNA', $mergingDependency, undef, $command, 'metrics/' , $workDirectory);
		$metricsJobId .= '$' .$rnaqcJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	
	##rawcount Matrix
	my $countDependency = undef;
	if($depends > 0) {
		$countDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'rawCounts'}}));
	}
	
	print "mkdir -p DGE\n";
	my $readCountDir = 'raw_counts' ;
	my $readcountExtension = '.readcounts.csv';
	my $outputDir = 'DGE';
	my $outputMatrix = 'rawCountMatrix.csv';
	$command = HtseqCount::refGtf2matrix($rH_cfg, LoadConfig::getParam($rH_cfg, 'htseq', 'referenceGtf'), $readCountDir, $readcountExtension, $outputDir, $outputMatrix);
	my $matrixJobId = undef;
	if(defined($command) && length($command) > 0) {
		$matrixJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'matrix', 'MATRIX', $countDependency, undef, $command, 'metrics/' , $workDirectory);
		$matrixJobId = '$' .$matrixJobId;
		$metricsJobId .= $matrixJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	
	##RPKM and Saturation
	my $countFile   = 'DGE/rawCountMatrix.csv';
	my $geneSizeFile     = LoadConfig::getParam($rH_cfg, 'saturation', 'geneSizeFile');
	my $rpkmDir = 'raw_counts';
	my $saturationDir = 'metrics';
	
	$command =  Metrics::saturation($rH_cfg, $countFile, $geneSizeFile, $rpkmDir, $saturationDir);
	my $saturationJobId = undef;
	if(defined($command) && length($command) > 0) {
		$saturationJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "saturation", undef, 'RPKM', $matrixJobId, undef, $command, 'metrics/' , $workDirectory);
		$metricsJobId .= '$' .$saturationJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	
	##fpkm Stats & Correlation
	my $fpkmDependency = undef;
	if($depends > 0) {
		$fpkmDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'fpkm'}}));
	}

	my $patern = '.fpkm_tracking';
	my $folder = 'fpkm/known';
	my $outputBaseName = 'metrics/fpkm';

	$command =  Metrics::fpkmCor($rH_cfg, $patern, $folder, $outputBaseName);
	my $fpkmJobId = undef;
	if(defined($command) && length($command) > 0) {
		$fpkmJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'fpkmStat', 'STATS_COR', $fpkmDependency, undef, $command, 'metrics/' , $workDirectory);
		$metricsJobId .= '$' .$fpkmJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	
	##readStats merge all files together and remove individuals ones
# 	my $mergeDependency = undef;
# 	if($depends > 0) {
# 		$mergeDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'),values(%{$globalDep{'aligning'}}));
# 	}
# 
# 	$patern = '.readstats.csv';
# 	$folder = 'metrics';
# 	$outputBaseName = 'metrics/readstats.AllSample.csv';
# 	
# 	$command =  Metrics::mergeReadStats($rH_cfg, $patern, $folder, $outputBaseName);
# 	my $mergeJobId = undef;
# 	if(defined($command) && length($command) > 0) {
# 		$mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'mergeRS', 'MERGEREADSTAT', $fpkmDependency, undef, $command, 'metrics/' , $workDirectory);
# 		$metricsJobId .= '$' .$mergeJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
# 	}
	$metricsJobId = substr $metricsJobId, 0, -1 ;
	return $metricsJobId;
}
 
sub dge {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'metrics'}{'metrics'};
	}
	
	print "mkdir -p DGE/output_jobs\n";
	
	my $countMatrix = 'DGE/rawCountMatrix.csv';
	my $outputDir = 'DGE';
	
	## edgeR
	my $command = DiffExpression::edgerPortable($rH_cfg, $designFilePath, $countMatrix, $outputDir);
	my $edgerJobId = undef;
	if(defined($command) && length($command) > 0) {
		$edgerJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "diffExpress", 'edger', 'EDGER', $jobDependency, undef, $command, 'DGE/' , $workDirectory);
		$edgerJobId = '$' .$edgerJobId;
	}
	
	## DESeq
	$command = DiffExpression::deseq($rH_cfg, $designFilePath, $countMatrix, $outputDir);
	my $deseqJobId = undef;
	if(defined($command) && length($command) > 0) {
		$deseqJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "diffExpress", 'deseq', 'DESEQ', $jobDependency, undef,$command, 'DGE/' , $workDirectory);
		$deseqJobId = '$' .$deseqJobId;
	}
	
	return $deseqJobId;
} 

sub goseq {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	
	my $cuffdiffDependency = undef;
	my $dgeDependency = undef;
	if($depends > 0) {
		$cuffdiffDependency = $globalDep{'cuffdiff'}{'cuffdiff'};
		$dgeDependency = $globalDep{'dge'}{'dge'};
	}

	my $columnsCuff = '2,13';
	my $columnsDge = '2,6';
	my $command;
	my $goseqJobId;
	for my $design (keys %{$rHoAoA_designGroup}) {
		## goseq for cuffdiff known results
		print "mkdir -p DGE/$design/output_jobs cuffdiff/$design/output_jobs\n";
		my $resultFileCuff = 'cuffdiff/known/' .$design .'/isoform_exp.diff' ;
		my $outputFileCuff = 'cuffdiff/known/' .$design .'/gene_ontology_results.csv';
		$command = DiffExpression::goseq($rH_cfg, $resultFileCuff, $outputFileCuff, $columnsCuff, '1');
		my $goCuffJobId = undef;
		if(defined($command) && length($command) > 0) {
			$goCuffJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "goseq","CUFFLINKS", 'GOCUFFDIFF' .$rH_jobIdPrefixe ->{$design} , $cuffdiffDependency, $design, $command, 'cuffdiff/' .$design, $workDirectory);
			$goseqJobId .= '$' .$goCuffJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
		## goseq for dge results
		my $resultFileDge = 'DGE/' .$design .'/dge_results.csv' ;
		my $outputFileDge = 'DGE/' .$design .'/gene_ontology_results.csv';
		$command = DiffExpression::goseq($rH_cfg, $resultFileDge, $outputFileDge, $columnsDge, '0');
		my $goDgeJobId = undef;
		if(defined($command) && length($command) > 0) {
			$goDgeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "goseq", "DESEQ", 'GODGE' .$rH_jobIdPrefixe ->{$design} , $dgeDependency, $design, $command, 'DGE/' .$design, $workDirectory);
			$goseqJobId .= '$' .$goDgeJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	}
	$goseqJobId = substr $goseqJobId, 0, -1 ;
	return $goseqJobId;
}

sub delivrable {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rHoAoA_designGroup   = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	
	my $goDependency ;
	if($depends > 0) {
		$goDependency = $globalDep{'goseq'}{'goseq'};
	}

}
1;
