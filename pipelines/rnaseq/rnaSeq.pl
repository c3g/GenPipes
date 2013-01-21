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
push(@steps, {'name' => 'metrics' , 'stepLoop' => 'sample' , 'output' => 'stats'});
push(@steps, {'name' => 'wiggle' , 'stepLoop' => 'sample' , 'output' => 'tracks'});
push(@steps, {'name' => 'rawCounts' , 'stepLoop' => 'sample' , 'output' => 'reads_count'});
push(@steps, {'name' => 'fpkm' , 'stepLoop' => 'sample' , 'output' => 'fpkm'});
push(@steps, {'name' => 'saturationRpkm' , 'stepLoop' => 'group' , 'output' => 'reads_count'});
push(@steps, {'name' => 'cuffdif' , 'stepLoop' => 'group' , 'output' => 'DGE'});
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
          my $jobIdVar = &$subref($current != ($opts{'s'}-1), \%cfg, undef, undef $rAoH_seqDictionary);
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
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $alignJobIdVarNameLane=undef;
    my $commands = Tophat::aln($rH_cfg, $sampleName, $rH_laneInfo, $rH_trimDetails->{'pair1'}, $rH_trimDetails->{'pair2'}, $rH_trimDetails->{'single1'}, $rH_trimDetails->{'single2'});
    if(defined $commands){
      my $alignJobIdVarNameLane = SubmitToCluster::printSubmitCmd($rH_cfg, "align", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'ALIGN', $jobDependency, $sampleName, $commands);
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

  my $command = Picard::merge($rH_cfg, $sampleName, $rAoH_sampleLanes);
  my $mergeJobId = undef;
  if(defined($command) && length($command) > 0) {
    $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "merging", undef, 'MERGELANES', $jobDependency, $sampleName, $command);
  }
  return $mergeJobId;
}

sub indelRealigner {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '$MERGELANES_JOB_ID';
  }

  print "mkdir -p $sampleName/realign\n";
  print "REALIGN_JOB_IDS=\"\"\n";
  my $processUnmapped = 1;
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $command = GATK::realign($rH_cfg, $sampleName, $seqName, $processUnmapped);
    my $intervalJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", $seqName, 'REALIGN', $jobDependency, $sampleName, $command);
    $intervalJobId = '$'.$intervalJobId;
    print 'REALIGN_JOB_IDS=${REALIGN_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$intervalJobId."\n";
    if($processUnmapped == 1) {
      $processUnmapped = 0;
    }
  }
  
  return '${REALIGN_JOB_IDS}';
}

sub mergeRealigned {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${REALIGN_JOB_IDS}';
  }

  my @seqNames;
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    push(@seqNames, $rH_seqInfo->{'name'});
  }
  my $command = Picard::mergeRealigned($rH_cfg, $sampleName, \@seqNames);
  my $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeRealigned", undef, 'MERGEREALIGN', $jobDependency, $sampleName, $command);
  return $mergeJobId;
}

sub fixmate {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MERGEREALIGN_JOB_ID}';
  }

  my $command = Picard::fixmate($rH_cfg, $sampleName);
  my $fixmateJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "fixmate", undef, 'FIXMATE', $jobDependency, $sampleName, $command);
  return $fixmateJobId;
}

sub markDup {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${FIXMATE_JOB_ID}';
  }

  my $command = Picard::markDup($rH_cfg, $sampleName);
  my $markDupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'MARKDUP', $jobDependency, $sampleName, $command);
  return $markDupJobId;
}

sub metrics {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MARKDUP_JOB_ID}';
  }

  my $command;
  my $bamFile = $sampleName.'/'.$sampleName.'.sorted.dup.bam';

  # Collect metrics
  $command = Picard::collectMetrics($rH_cfg, $sampleName);
  my $collectMetricsJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "collectMetrics", undef, 'COLLECTMETRICS', $jobDependency, $sampleName, $command);
  
  # Compute genome coverage
  $command = GATK::genomeCoverage($rH_cfg, $sampleName);
  my $genomeCoverageJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "genomeCoverage", undef, 'GENOMECOVERAGE', $jobDependency, $sampleName, $command);

  # Compute CCDS coverage
  my $outputPrefix = $sampleName.'/'.$sampleName.'.sorted.dup.CCDS.coverage';
  $command = GATK::targetCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  my $targetCoverageJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "targetCoverage", undef, 'TARGETCOVERAGE', $jobDependency, $sampleName, $command);

  # Generate IGV track
  $command = IGVTools::computeTDF($rH_cfg, $sampleName, $bamFile);
  my $igvtoolsTDFJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "computeTDF", undef, 'IGVTOOLS', $jobDependency, $sampleName, $command);

  # Compute flags
  my $output = $sampleName.'/'.$sampleName.'.sorted.dup.bam.flagstat';
  $command = SAMTools::flagstat($rH_cfg, $sampleName, $bamFile, $output);
  my $flagstatJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "flagstat", undef, 'FLAGSTAT', $jobDependency, $sampleName, $command);

  # return the longest one...not ideal
  return $genomeCoverageJobId;
}

sub sortQname {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MARKDUP_JOB_ID}';
  }
}

#push(@steps, {'name' => 'countTelomere'});
sub fullPileup {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MARKDUP_JOB_ID}';
  }

  my $bamFile = $sampleName.'/'.$sampleName.'.sorted.dup.bam';
  print 'mkdir -p '.$sampleName.'/mpileup/'."\n";
  print "RAW_MPILEUP_JOB_IDS=\"\"\n";
  my $catCommand = 'zcat ';
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $outputPerSeq = $sampleName.'/mpileup/'.$sampleName.'.'.$seqName.'.mpileup.gz';
    my $command = SAMtools::rawmpileup($rH_cfg, $sampleName, $bamFile, $seqName, $outputPerSeq);
    my $mpileupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "rawmpileup", $seqName, 'RAW_MPILEUP', $jobDependency, $sampleName, $command);
    $mpileupJobId = '$'.$mpileupJobId;
    print 'RAW_MPILEUP_JOB_IDS=${RAW_MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$mpileupJobId."\n";

    $catCommand .= $outputPerSeq .' ';
  }

  my $output = $sampleName.'/mpileup/'.$sampleName.'.mpileup.gz';
  $catCommand .= '| gzip -c --best > '.$output;

  my $catJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "rawmpileup_cat", undef, 'RAW_MPILEUP_CAT', undef, "\$RAW_MPILEUP_JOB_IDS", $catCommand);
  return $catJobId;
}

1;
