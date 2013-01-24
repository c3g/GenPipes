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


use GATK;
use IGVTools;
use LoadConfig;
use Picard;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SubmitToCluster;
use TophatBowtie;
use Trimmomatic;
#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'trimming' , 'stepLoop' => 'sample' , 'output' => 'reads'});
push(@steps, {'name' => 'aligning' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
push(@steps, {'name' => 'merging' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
push(@steps, {'name' => 'wiggle' , 'stepLoop' => 'sample' , 'output' => 'tracks'});
push(@steps, {'name' => 'rawCounts' , 'stepLoop' => 'sample' , 'output' => 'reads_count'});
push(@steps, {'name' => 'fpkm' , 'stepLoop' => 'sample' , 'output' => 'fpkm'});
push(@steps, {'name' => 'metrics' , 'stepLoop' => 'group' , 'output' => 'stats'});
#push(@steps, {'name' => 'saturationRpkm' , 'stepLoop' => 'group' , 'output' => 'reads_count'}); included in metrics
push(@steps, {'name' => 'cuffdiff' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'dge' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'delivrable' , 'stepLoop' => 'group' , 'output' => 'Delivrable'});


my %globalDep;
for my $stepName (@steps) { 
	$globalDep{$stepName -> {'name'} } ={};
}



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
  
  if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'e'}) || !defined($opts{'n'}) || !defined($opts{'d'}|| !defined($opts{'w'} ) ) {
    printUsage();
    exit(1);
  }

  my %cfg = LoadConfig->readConfigFile($opts{'c'});
  my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($opts{'n'});
  my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);
  my $designFilePath = $opts{'d'};
  my $workDirectory = $opts{'w'};

  my $latestBam;

    
    for(my $current = $opts{'s'}-1; $current <= ($opts{'e'}-1); $current++) {
       my $fname = $steps[$current]->{'name'};
       my $loopType = $steps[$current]->{'stepLoop'};
       my $outputStep = $steps[$current]->{'output'};
       my $subref = \&$fname;
       if ($loopType == 'sample') {
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
          my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
	  my $outputLocation = $outputStep. "/" .$sampleName
	  SubmitToCluster::initSubmit(\%cfg, $outputLocation);
          # Tests for the first step in the list. Used for dependencies.
          my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, $sampleName, $rAoH_sampleLanes, $rAoH_seqDictionary);
	  $globalDep{$fname}{$sampleName -> {$jobIdVar}};
        }
       }
       else {
	SubmitToCluster::initSubmit(\%cfg, $outputLocation);
          # Tests for the first step in the list. Used for dependencies.
          my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, $rHoAoH_sampleInfo, undef $rAoH_seqDictionary);
	  $globalDep{$fname}{$fname -> {$jobIdVar}};
       }
    }  
}

sub trimming {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $trimJobIdVarNameSample = undef;
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo);
    my $trimJobIdVarNameLane=undef;
    if(length($rH_trimDetails->{'command'}) > 0) {
      $trimJobIdVarNameLane = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM', undef, $sampleName, $rH_trimDetails->{'command'}, $workDirectory);
      $trimJobIdVarNameSample .= '$'.$trimJobIdVarNameLane .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
      #TODO calcReadCounts.sh
    }
   }
    $trimJobIdVarNameSample = substr $trimJobIdVarNameSample 0 -1;
  return $trimJobIdVarNameSample;	
}

sub aligning {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $alignJobIdVarNameSample = undef;
  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = $globalDep{'trimming'}{$sampleName};
  }
  
  print "mkdir -p alignment\n";
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $alignJobIdVarNameLane=undef;
    my $commands = Tophat::aln($rH_cfg, $sampleName, $rH_laneInfo, $rH_trimDetails->{'pair1'}, $rH_trimDetails->{'pair2'}, $rH_trimDetails->{'single1'}, $rH_trimDetails->{'single2'});
    if(defined $commands && length($command) > 0){
      my $alignJobIdVarNameLane = SubmitToCluster::printSubmitCmd($rH_cfg, "align", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'ALIGN', $jobDependency, $sampleName, $commands, $workDirectory);
      $alignJobIdVarNameSample .= '$'. $alignJobIdVarNameLane .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'); 
    } 
    $alignJobIdVarNameSample = substr $alignJobIdVarNameSample 0 -1;
  }
  return $alignJobIdVarNameSample;
}

sub merging {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = $globalDep{'aligning'}{$sampleName};
  }

  ##Merging
  my $inputBAM ; 
  my $outputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.bam" ;
  my @alignFiles;
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $laneDirectory = "alignment/" . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    $inputBAM = $laneDirectory . 'accepted_hits.bam';
    push(@alignFiles, $inputBAM) ;
  }
  my $command = Picard::mergeFiles($rH_cfg, $sampleName, $rAoH_sampleLanes, $outputBAM);
  my $mergeJobId = undef;
  if(defined($command) && length($command) > 0) {
    $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergFiles", undef, 'MERGELANES', $jobDependency, $sampleName, $command);
    $mergeJobId = '$'.$mergeJobId;
  }
  ## reorder
  $inputBAM = $outputBAM
  $outputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.karyotypic.bam";
  $command = Picard::mergeFiles($rH_cfg, $sampleName, $inputBAM, $outputBAM);
  my $reorderJobId = undef;
  if(defined($command) && length($command) > 0) {
    $reorderJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "reorderSam", undef, 'REORDER', $mergeJobId, $sampleName, $command);
    $reorderJobId = '$'.$reorderJobId;
  }
  ## mark duplicates
  $inputBAM = $outputBAM
  $outputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.bam"
  my $duplicatesMetricsFile = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.metrics";
  $command = Picard::markDup($rH_cfg, $sampleName, $inputBAM, $outputBAM);
  my $markDupJobId = undef;
  if(defined($command) && length($command) > 0) {
    $markDupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'MARKDUP', $reorderJobId, $sampleName, $command);
    $markDupJobId = '$'.$markDupJobId
  }
  return $markDupJobId;
}

sub wiggle {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}

	my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.bam" ; 
	#testing for strand-specificity

	my $strandSPecificityInfo = LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
	my @strandJobId ;
	if($strandSPecificityInfo != "fr-unstranded") {
	## strand specific 
		my @outputBAM = {'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.forward.bam' ,  'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.reverse.bam'};
		my @command = Wiggle::strandBam($rH_cfg, $sampleName, $inputBAM, \@outputBAM);
		if(defined(@command) && length(@command) > 1) { 
			my $strandJobIdF = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", 'FORWARD', 'STRANDSPEC', $jobDependency, $sampleName, @command->[0]);
			push(@strandJobId, '$'.$strandJobIdF );
			my $strandJobIdR = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", 'REVERSE', 'STRANDSPEC', $jobDependency, $sampleName, @command->[1]);
			push(@strandJobId, '$'.$strandJobIdR );
		}
		my @outputBedGraph = {'alignment/' . $sampleName . '/' . $sampleName . '.forward.bedGraph' ,  'alignment/' . $sampleName . '/' . $sampleName . '.reverse.bedGraph'};
		my @outputWiggle = {'alignment/' . $sampleName . '/' . $sampleName . '.forward.bw' ,  'alignment/' . $sampleName . '/' . $sampleName . '.reverse.bw'};
		my @prefixJobName = { 'FORWARD', 'REVERSE'};
	}
	else {
		my @outputBAM = {$inputBAM};
		push(@strandJobId, $jobDependency);
		my @outputBedGraph = {'alignment/' . $sampleName . '/' . $sampleName . '.bedGraph'};
		my @outputWiggle = {'alignment/' . $sampleName . '/' . $sampleName . '.bw' };
		my @prefixJobName = { undef } ;
	}
	my $wiggleJobId ;
	for(my $i = 0; $i <@outputBAM; $i++) {
		my $command = Wiggle::strandBam($rH_cfg, $sampleName, $inputBAM, @outputBedGraph->[$i], @outputWiggle->[$i]);
		if(defined($command) && length($command) > 0) {
			my $tmpJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", @prefixJobName->[$i], 'WIGGLE', @strandJobId->[$i], $sampleName, $command);
			$wiggleJobId .= '$'.$tmpJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
		}
	} 
	$wiggleJobId = substr $wiggleJobId 0 -1;
	return $wiggleJobId;	
}

sub rawCounts {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}
	print "mkdir -p raw_counts\n";
	my $inputBam = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam' ;
        my $sortedBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.queryNameSorted.bam' ;
	my $inputGtf = LoadConfig::getParam($rH_cfg, 'htseq', 'referenceGtf');
	my $outputCount = 'raw_counts/' . $sampleName . '.readcounts.csv';
	my $sortOrder = 'queryname';
	my $strandInfo;
	my $strandSPecificityInfo = LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
	if($strandSPecificityInfo != "fr-unstranded") {
		 $strandInfo= 'yes';
	}
	else {
		$strandInfo= 'no';
	}
	## query sort the bam
	my $sortJobId;
	my $command = Picard::sortSam($rH_cfg, $sampleName, $inputBAM, $sortedBAM, $sortOrder);
	if(defined($command) && length($command) > 0) {
		$sortJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "sortSam", undef, 'QNSORT', $jobDependency, $sampleName, $command);
		$sortJobId='$'.$sortJobId
	}
	## count reads
        my $countJobId;
	my $command = HtseqCount::readCount($rH_cfg, $sortedBAM, $inputGtf, $outputCount, $strandInfo); 
	if(defined($command) && length($command) > 0) {
		$countJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "htseq", undef, 'RAWCOUNT', $sortJobId, $sampleName, $command);
		$countJobId='$'.$countJobId
	}
	return $countJobId;
}

sub fpkm {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = $globalDep{'merging'}{$sampleName};
	}
	
	print "mkdir -p fpkm/known fpkm/denovo\n";
	my $inputBam = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam' ;
	my $outputKnown = 'fpkm/known/' . $sampleName;
	my $outputDeNovo = 'fpkm/denovo' . $sampleName;
	my $gtfOption = '-G ' .LoadConfig::getParam($rH_cfg, 'fpkm','referenceGtf');
	
	
	## known FPKM
	my $fpkmJobId;
	my $command = Cufflinks::fpkm($rH_cfg, $sampleName, $inputBAM, $outputKnown, $gtfOption);
	if(defined($command) && length($command) > 0) {
		my $fpkmKnownJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "fpkm", "KNOWN", 'FPKM', $jobDependency, $sampleName, $command);
		$fpkmJobId='$'.$fpkmKnownJobId .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
	}
	## denovo FPKM
	my $command = Cufflinks::fpkm($rH_cfg, $sampleName, $inputBAM, $outputDeNovo, undef);
	if(defined($command) && length($command) > 0) {
		my $fpkmDeNovoJobId=SubmitToCluster::printSubmitCmd($rH_cfg, "fpkm", "DENOVO", 'FPKM', $jobDependency, $sampleName, $command);
		$fpkmJobId .='$' .$fpkmDeNovoJobId;
	}
	return $fpkmJobId;
}

sub cuffdiff {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	
	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = join(':',values(%{$globalDep ->{'fpkm'}}));
	}
	print "mkdir -p cuffdiff/Known cuffdiff/denovo\n";
	##get design groups
	my $rHoAoA_command = Cufflinks::getDesign($rH_cfg,$designFilePath);
}

sub metrics {
	my $depends = shift;
	my $rH_cfg = shift;
	my $rHoAoH_sampleInfo = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;

	my $jobDependency = undef;
	if($depends > 0) {
		$jobDependency = join(':',values(%{$globalDep ->{'merging'}}));
	}
	## RNAseQC metrics
	print "echo -e \"Sample\tBamFile\tNote\" >  alignment/rnaseqc.samples.txt\n"
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
		print 'echo -e \"' .$sampleName .'\talignment/' .$sampleName. '/'. $sampleName .'.merged.mdup.bam\t' .LoadConfig::getParam($rH_cfg, 'metricsRNA', 'projectName'). '\" >>  alignment/rnaseqc.samples.txt'
	}
	
	print "mkdir -p metrics\n";
	my $sampleList = "alignment/rnaseqc.samples.txt"
	my $outputFolder = "metrics"
	my $command = Metrics::rnaQc($rH_cfg, $sampleList, $outputFolder);
	my $metricsJobId = undef;
	if(defined($command) && length($command) > 0) {
		$metricsJobId= SubmitToCluster::printSubmitCmd($rH_cfg, "metricsRNA", undef, 'METRICSRNA', $jobDependency, undef, $command);
	}
	


	##Saturation


	##fpkm Stats & Correlation


	##BamContent

	return $metrics;
}
  

1;
