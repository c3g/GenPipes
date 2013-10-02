#!/usr/bin/perl

=head1 NAME

I<dnaSeq>

=head1 SYNOPSIS

dnaSeq.pl

=head1 DESCRIPTION

B<dnaSeq> Is the main variant discovery pipeline.

=head1 AUTHOR

B<Louis Letourneau> - I<louis.letourneau@mail.mcgill.ca>

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
use POSIX;

use BWA;
use GATK;
use IGVTools;
use LoadConfig;
use Picard;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SnpEff;
use SubmitToCluster;
use Trimmomatic;
use VCFtools;
use Metrics;
#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'trimAndAlign', 'stepLoop' => 'sample', 'parentStep' => undef});
push(@steps, {'name' => 'metricsLanes', 'stepLoop' => 'experiment', 'parentStep' => 'trimAndAlign'});
push(@steps, {'name' => 'mergeLanes', 'stepLoop' => 'sample', 'parentStep' => 'trimAndAlign'});
push(@steps, {'name' => 'indelRealigner', 'stepLoop' => 'sample', 'parentStep' => 'mergeLanes'});
push(@steps, {'name' => 'mergeRealigned', 'stepLoop' => 'sample', 'parentStep' => 'indelRealigner'});
push(@steps, {'name' => 'fixmate', 'stepLoop' => 'sample', 'parentStep' => 'mergeRealigned'});
push(@steps, {'name' => 'markDup', 'stepLoop' => 'sample', 'parentStep' => 'fixmate'});
push(@steps, {'name' => 'recalibration', 'stepLoop' => 'sample', 'parentStep' => 'markDup'});
push(@steps, {'name' => 'metrics', 'stepLoop' => 'sample', 'parentStep' => 'recalibration'});
push(@steps, {'name' => 'metricsLibrarySample', 'stepLoop' => 'experiment', 'parentStep' => 'metrics'});
push(@steps, {'name' => 'fullPileup', 'stepLoop' => 'sample', 'parentStep' => 'recalibration'});
push(@steps, {'name' => 'snpAndIndelBCF', 'stepLoop' => 'experiment', 'parentStep' => 'recalibration'});
push(@steps, {'name' => 'mergeFilterBCF', 'stepLoop' => 'experiment', 'parentStep' => 'snpAndIndelBCF'});
#push(@steps, {'name' => 'filterNStretches', 'stepLoop' => 'experiment', 'parentStep' => 'mergeFilterBCF'});
push(@steps, {'name' => 'flagMappability', 'stepLoop' => 'experiment', 'parentStep' => 'mergeFilterBCF'});
push(@steps, {'name' => 'snpIDAnnotation', 'stepLoop' => 'experiment', 'parentStep' => 'flagMappability'});
#push(@steps, {'name' => 'snpEffect', 'stepLoop' => 'experiment', 'parentStep' => 'snpIDAnnotation'});
push(@steps, {'name' => 'metricsSNV', 'stepLoop' => 'experiment', 'parentStep' => 'snpIDAnnotation'});
# push(@steps, {'name' => 'breakDancer', 'stepLoop' => 'experiment', 'parentStep' => 'recalibration'});
# push(@steps, {'name' => 'pindel', 'stepLoop' => 'experiment', 'parentStep' => 'breakDancer'});
# push(@steps, {'name' => 'cnv', 'stepLoop' => 'experiment', 'parentStep' => 'recalibration'});
# push(@steps, {'name' => 'metricsSV', 'stepLoop' => 'experiment', 'parentStep' => ('cnv','pindel')});
push(@steps, {'name' => 'deliverable' , 'stepLoop' => 'experiment' , 'parentStep' => ('metricsLanes','metricsSample','metricsSNV')});


my %globalDep;
for my $stepName (@steps) {
  $globalDep{$stepName -> {'name'} } ={};
}

&main();

sub printUsage {
  print "\nUsage: perl ".$0." -c config.ini -s start -e end -n SampleSheet.csv\n";
  print "\t-c  config file\n";
  print "\t-s  start step, inclusive\n";
  print "\t-e  end step, inclusive\n";
  print "\t-n  nanuq sample sheet\n";
  print "\n";
  print "Steps:\n";
  for(my $idx=0; $idx < @steps; $idx++) {
    print "".($idx+1).'- '.$steps[$idx]->{'name'}."\n";
  }
  print "\n";
}

sub main {
  my %opts;
  getopts('c:s:e:n:', \%opts);
  
  if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'e'}) || !defined($opts{'n'})) {
    printUsage();
    exit(1);
  }

  my %cfg = LoadConfig->readConfigFile($opts{'c'});
  my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($opts{'n'});
  my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);

  my $latestBam;
  my @sampleNames = keys %{$rHoAoH_sampleInfo};

  print STDERR "Samples: ".scalar(@sampleNames)."\n";

  SubmitToCluster::initPipeline;

  my $currentStep;
  for(my $idx=0; $idx < @sampleNames; $idx++){
    my $sampleName = $sampleNames[$idx];
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};

    for($currentStep = $opts{'s'}-1; $currentStep <= ($opts{'e'}-1); $currentStep++) {
      my $fname = $steps[$currentStep]->{'name'};
      my $subref = \&$fname;

      if ($steps[$currentStep]->{'stepLoop'} eq 'sample') {
        # Tests for the first step in the list. Used for dependencies.
        my $jobIdVar = &$subref($currentStep, \%cfg, $sampleName, $rAoH_sampleLanes, $rAoH_seqDictionary); 
        $globalDep{$fname}->{$sampleName} = $jobIdVar;
      }
    }
  }  

  for($currentStep = $opts{'s'}-1; $currentStep <= ($opts{'e'}-1); $currentStep++) {
    if($steps[$currentStep]->{'stepLoop'} eq 'experiment') {
      my $fname = $steps[$currentStep]->{'name'};
      my $subref = \&$fname;

      my $jobIdVar = &$subref($currentStep, \%cfg, $rHoAoH_sampleInfo, $rAoH_seqDictionary);
      $globalDep{$fname}->{'experiment'} = $jobIdVar;
    }
  }
}

sub trimAndAlign {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  print "BWA_JOB_IDS=\"\"\n";
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $rgId = $rH_laneInfo->{'libraryBarcode'} . "_" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rgSampleName = $rH_laneInfo->{'name'};
    my $rgLibrary = $rH_laneInfo->{'libraryBarcode'};
    my $rgPlatformUnit = $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rgCenter = LoadConfig::getParam( $rH_cfg, 'aln', 'bwaInstitution' );

    my $outputDir = 'reads/'.$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    print 'mkdir -p '.$outputDir."\n";
    my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo, $outputDir);
    my $trimJobIdVarName=undef;
    if(length($rH_trimDetails->{'command'}) > 0) {
      $trimJobIdVarName = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM', undef, $sampleName, $rH_trimDetails->{'command'});
      $trimJobIdVarName = '$'.$trimJobIdVarName;
    }

    my $outputAlnDir = 'alignment/'.$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    print 'mkdir -p '.$outputAlnDir."\n";
    my $outputAlnPrefix = $outputAlnDir.'/'.$sampleName;
    my $rA_commands = BWA::aln($rH_cfg, $sampleName, $rH_trimDetails->{'pair1'}, $rH_trimDetails->{'pair2'}, $rH_trimDetails->{'single1'}, $outputAlnPrefix, $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter);
    if(defined($rA_commands)) {
      if(@{$rA_commands} == 3) {
        my $read1JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read1.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ1ALN', $trimJobIdVarName, $sampleName, $rA_commands->[0]);
        $read1JobId = '$'.$read1JobId;
        my $read2JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read2.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ2ALN', $trimJobIdVarName, $sampleName, $rA_commands->[1]);
        $read2JobId = '$'.$read2JobId;
        my $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'sampe.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA', $read1JobId.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$read2JobId, $sampleName, $rA_commands->[2]);
        $bwaJobId = '$'.$bwaJobId;
        print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$bwaJobId."\n";
      }
      else {
        my $readJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READALN', $trimJobIdVarName, $sampleName, $rA_commands->[0]);
        $readJobId = '$'.$readJobId;
        my $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'samse.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA',  $readJobId, $sampleName, $rA_commands->[1]);
        $bwaJobId = '$'.$bwaJobId;
        print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$bwaJobId."\n";
      } 
    }
  }
  return '$BWA_JOB_IDS';
}


sub metricsLanes {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;
  
  my $libraryType =  undef;
  my $fkey =  (keys %{$rHoAoH_sampleInfo})[0] ;
  my @fvals = @{$rHoAoH_sampleInfo->{$fkey}};
  my $finfo = $fvals[0];
  if ( $finfo->{'runType'} eq "SINGLE_END" ) {
    $libraryType = 'single';
  } elsif ($finfo->{'runType'} eq "PAIRED_END" ) {
    $libraryType = 'paired';
  }
  my $trimmingDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};

  my @sampleNames = keys %{$rHoAoH_sampleInfo};
  my $jobDependencies = "";
  for(my $idx=0; $idx < @sampleNames; $idx++){
    my $sampleName = $sampleNames[$idx];
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
    if(defined($globalDep{$parentStep}->{$sampleName})){
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$globalDep{$parentStep}->{$sampleName};
    }
  }
  if(length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }
  $trimmingDependency = $jobDependencies;
  
  my $folder = 'reads';
  my $pattern = 'trim.stats.csv';
  my $ouputFile = 'metrics/trimming.stats';
  print "mkdir -p metrics\n";
  my $command;
  $command = Metrics::mergeTrimmomaticStats($rH_cfg,  $libraryType, $pattern, $folder, $ouputFile);
  my $metricsJobId = undef;
  if(defined($command) && length($command) > 0) {
          my $trimMetricsJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "trimMetrics", undef, 'TRIMMETRICS', $trimmingDependency, undef, $command, LoadConfig::getParam($rH_cfg, "default", 'sampleOutputRoot') );
          $metricsJobId = '$' .$trimMetricsJobId;
  }
  return $metricsJobId;
}


sub mergeLanes {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

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
  if(defined($command) && length($command) > 0) {
    $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeLanes", undef, 'MERGELANES', $jobDependency, $sampleName, $command);
  }
  else {
    return undef;
  }
  return '$'.$mergeJobId;
}

sub indelRealigner {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }
  print 'mkdir -p alignment/'.$sampleName."/realign\n";

  my $nbRealignJobs = LoadConfig::getParam( $rH_cfg, 'indelRealigner', 'nbRealignJobs' );
  if($nbRealignJobs > 50) {
    warn "Number of realign jobs is >50. This is usually much. Anything beyond 20 can be problematic.\n";
  }

  if($nbRealignJobs <= 1) {
    my $command = GATK::realign($rH_cfg, $sampleName, 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam', undef, 'alignment/'.$sampleName.'/realign/all');
    my $intervalJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", undef, 'REALIGN', $jobDependency, $sampleName, $command);
    print 'REALIGN_JOB_IDS=$'.$intervalJobId."\n";
  }
  else {
    #Keep space for the exclude realignment at the end.
    $nbRealignJobs--;
    my @chrToProcess;
    for (my $idx=0; $idx < $nbRealignJobs; $idx++) {
      push(@chrToProcess, $rAoH_seqDictionary->[$idx]->{'name'});
    }

    print "REALIGN_JOB_IDS=\"\"\n";
    my $processUnmapped = 1;
    my @excludeList;
    for my $seqName (@chrToProcess) {
      push(@excludeList, $seqName);
      my $command = GATK::realign($rH_cfg, $sampleName, 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam', $seqName, 'alignment/'.$sampleName.'/realign/'.$seqName, $processUnmapped);
      my $intervalJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", $seqName, 'REALIGN', $jobDependency, $sampleName, $command);
      $intervalJobId = '$'.$intervalJobId;
      if($processUnmapped == 1) {
        print 'REALIGN_JOB_IDS='.$intervalJobId."\n";
        $processUnmapped = 0;
      }
      else {
        print 'REALIGN_JOB_IDS=${REALIGN_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$intervalJobId."\n";
      }
    }

    my $command = GATK::realign($rH_cfg, $sampleName, 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam', undef, 'alignment/'.$sampleName.'/realign/others', 1, \@excludeList);
    my $intervalJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", 'others', 'REALIGN', $jobDependency, $sampleName, $command);
    print 'REALIGN_JOB_IDS=${REALIGN_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').'$'.$intervalJobId."\n";
  }
  return '${REALIGN_JOB_IDS}';
}

sub mergeRealigned {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $latestBam;
  my @inputBams;
  my $outputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.realigned.qsorted.bam';

  my $command;
  my $nbRealignJobs = LoadConfig::getParam( $rH_cfg, 'indelRealigner', 'nbRealignJobs' );
  if($nbRealignJobs <= 1) {
    my $command = 'echo "ln -s alignment/'.$sampleName.'/realign/all.bam '.$outputBAM.'"';
  }
  else {
    #Keep space for the exclude realignment at the end.
    $nbRealignJobs--;
    my @chrToProcess;
    for (my $idx=0; $idx < $nbRealignJobs; $idx++) {
      my $input = 'alignment/'.$sampleName.'/realign/'.$rAoH_seqDictionary->[$idx]->{'name'}.'.bam';
      push(@inputBams, $input);
    }
    push(@inputBams, 'alignment/'.$sampleName.'/realign/others.bam');

    $command = Picard::mergeFiles($rH_cfg, $sampleName, \@inputBams, $outputBAM);
  }

  my $mergeJobId = undef;
  if(defined($command) && length($command) > 0) {
    $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeRealign", undef, 'MERGEREALIGN', $jobDependency, $sampleName, $command);
  }

  if(defined($mergeJobId)){
    return '$'.$mergeJobId;
  }
  return undef;
}

sub fixmate {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $inputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.realigned.qsorted.bam';
  my $outputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.matefixed.sorted.bam';

  my $command = Picard::fixmate($rH_cfg, $inputBAM, $outputBAM);
  my $fixmateJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "fixmate", undef, 'FIXMATE', $jobDependency, $sampleName, $command);
  return '$'.$fixmateJobId;
}

sub markDup {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $inputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.matefixed.sorted.bam';
  my $outputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.bam';
  my $outputMetrics = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.metrics';

  my $command = Picard::markDup($rH_cfg, $sampleName, $inputBAM, $outputBAM, $outputMetrics);
  my $markDupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'MARKDUP', $jobDependency, $sampleName, $command);
  return '$'.$markDupJobId;
}

sub recalibration {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $inputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.bam';
  my $outputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup';

  my $command = GATK::recalibration($rH_cfg, $sampleName, $inputBAM, $outputBAM);
  my $recalibrationJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "recalibration", undef, 'RECAL', $jobDependency, $sampleName, $command);
  return '$'.$recalibrationJobId;
}

sub metrics {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $command;
  my $bamFile = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.bam';

  # Collect metrics
  my $outputMetrics = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.all.metrics';
  $command = Picard::collectMetrics($rH_cfg, $bamFile, $outputMetrics);
  my $collectMetricsJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "collectMetrics", undef, 'COLLECTMETRICS', $jobDependency, $sampleName, $command);
  
  # Compute genome coverage
  my $outputPrefix = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.all.coverage';
  $command = GATK::genomeCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  my $genomeCoverageJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "genomeCoverage", undef, 'GENOMECOVERAGE', $jobDependency, $sampleName, $command);

  # Compute CCDS coverage
  $outputPrefix = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.CCDS.coverage';
  $command = GATK::targetCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  my $targetCoverageJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "targetCoverage", undef, 'TARGETCOVERAGE', $jobDependency, $sampleName, $command);

  # Generate IGV track
  $command = IGVTools::computeTDF($rH_cfg, $bamFile);
  my $igvtoolsTDFJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "computeTDF", undef, 'IGVTOOLS', $jobDependency, $sampleName, $command);

  # Compute flags
  my $output = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.bam.flagstat';
  $command = SAMtools::flagstat($rH_cfg, $bamFile, $output);
  my $flagstatJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "flagstat", undef, 'FLAGSTAT', $jobDependency, $sampleName, $command);

  # return the longest one...not ideal
  return '$'.$genomeCoverageJobId;
}


sub metricsLibrarySample {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;
  
  
  my $metricsDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};

  my @sampleNames = keys %{$rHoAoH_sampleInfo};
  my $jobDependencies = "";
  for(my $idx=0; $idx < @sampleNames; $idx++){
    my $sampleName = $sampleNames[$idx];
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
    if(defined($globalDep{$parentStep}->{$sampleName})){
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$globalDep{$parentStep}->{$sampleName};
    }
  }
  if(length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }
  $metricsDependency = $jobDependencies;
  
  my $folder = 'alignment/';
  my $ouputFile = 'metrics/SampleMetrics.stats';
  my $experimentType = LoadConfig::getParam($rH_cfg, 'default', 'experimentType');
  print "mkdir -p metrics\n";
  my $command;
  $command = Metrics::mergeSampleDnaStats($rH_cfg,  $experimentType, $folder, $ouputFile);
  my $metricsJobId = undef;
  if(defined($command) && length($command) > 0) {
          my $trimMetricsJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "sampleMetrics", undef, 'SAMPLEMETRICS', $metricsDependency, undef, $command, LoadConfig::getParam($rH_cfg, "default", 'sampleOutputRoot') );
          $metricsJobId = '$' .$trimMetricsJobId;
  }
  return $metricsJobId;
}


sub sortQname {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }
}

#push(@steps, {'name' => 'countTelomere'});
sub fullPileup {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $bamFile = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.bam';
  my $outputDir = 'alignment/'.$sampleName.'/mpileup/';

  print 'mkdir -p '.$outputDir."\n";
  print "RAW_MPILEUP_JOB_IDS=\"\"\n";
  my $catCommand = 'zcat ';
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $outputPerSeq = $outputDir.$sampleName.'.'.$seqName.'.mpileup.gz';
    my $command = SAMtools::rawmpileup($rH_cfg, $sampleName, $bamFile, $seqName, $outputPerSeq);
    my $mpileupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "rawmpileup", $seqName, 'RAW_MPILEUP', $jobDependency, $sampleName, $command);
    $mpileupJobId = '$'.$mpileupJobId;
    print 'RAW_MPILEUP_JOB_IDS=${RAW_MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$mpileupJobId."\n";

    $catCommand .= $outputPerSeq .' ';
  }

  my $output = $outputDir.$sampleName.'.mpileup.gz';
  $catCommand .= '| gzip -c --best > '.$output;

  my $catJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "rawmpileup_cat", undef, 'RAW_MPILEUP_CAT', "\$RAW_MPILEUP_JOB_IDS", $sampleName, $catCommand);
  return '$'.$catJobId;
}

sub snpAndIndelBCF {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $parentStep = $steps[$stepId]->{'parentStep'};

  my @sampleNames = keys %{$rHoAoH_sampleInfo};
  my $jobDependencies = "";
  my @inputFiles;
  for(my $idx=0; $idx < @sampleNames; $idx++){
    my $sampleName = $sampleNames[$idx];
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
    if(defined($globalDep{$parentStep}->{$sampleName})){
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$globalDep{$parentStep}->{$sampleName};
    }
    push(@inputFiles, 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.bam');
  }
  if(length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }

  my $outputDir = 'variants/rawBCF/';
  print 'mkdir -p '.$outputDir."\n";
  print "MPILEUP_JOB_IDS=\"\"\n";

  my $nbJobs = LoadConfig::getParam( $rH_cfg, 'mpileup', 'approxNbJobs' );
  my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);
  my $isFirst=1;
  for my $region (@{$rA_regions}) {
    my $command = SAMtools::mpileup($rH_cfg, 'allSamples', \@inputFiles, $region, $outputDir);
    $region =~ s/:/_/g;
    my $mpileupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $region, 'MPILEUP', $jobDependencies, 'allSamples', $command);
    $mpileupJobId = '$'.$mpileupJobId;
    if($isFirst==1) {
      print 'MPILEUP_JOB_IDS='.$mpileupJobId."\n";
      $isFirst=0;
    }
    else {
      print 'MPILEUP_JOB_IDS=${MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$mpileupJobId."\n";
    }
  }

  return '${MPILEUP_JOB_IDS}';
}

sub generateApproximateWindows {
  my $nbJobs = shift;
  my $rAoH_seqDictionary = shift;

  my @retVal;
  if($nbJobs <= scalar(@{$rAoH_seqDictionary})) {
    for my $rH_seqInfo (@$rAoH_seqDictionary) {
      push(@retVal, $rH_seqInfo->{'name'}.':1-'.$rH_seqInfo->{'size'});
    }
  }
  else{
    $nbJobs -= @$rAoH_seqDictionary;
    my $totalSize = 0;
    for my $rH_seqInfo (@$rAoH_seqDictionary) {
      $totalSize += $rH_seqInfo->{'size'};
    }
    my $approxWindow = floor($totalSize / $nbJobs);

    for my $rH_seqInfo (@$rAoH_seqDictionary) {
      for(my $idx=1; $idx <= $rH_seqInfo->{'size'}; $idx += $approxWindow) {
        my $end = $idx+$approxWindow-1;
        if($end > $rH_seqInfo->{'size'}) {
          $end = $rH_seqInfo->{'size'};
        }

        my $region = $rH_seqInfo->{'name'}.':'.$idx.'-'.$end;
        push(@retVal, $region);
      }
    }
  }

#  my @retVal;
#  $totalSize = 0;
#  my $currentregion = "";
#  my $approxWindowRemaining = $approxWindow;
#  for my $rH_seqInfo (@$rAoH_seqDictionary) {
#    for(my $idx=1; $idx <= $rH_seqInfo->{'size'}; $idx += $approxWindowRemaining) {
#      my $end = $idx+$snvWindow-1;
#      my $hitEnd = 0;
#      if($end > $rH_seqInfo->{'size'}) {
#        $end = $rH_seqInfo->{'size'};
#        $hitEnd = 1;
#        $approxWindowRemaining -= ($end - $idx) +1;
#      }
#
#      my $region = $rH_seqInfo->{'name'}.':'.$idx.'-'.$end;
#      if(length($currentregion) == 0) {
#        $currentregion = $region;
#      }
#      else {
#        $currentregion .= ','.$region;
#      }
#
#      if($hitEnd == 0) {
#        push(@retVal, $currentregion);
#        $currentregion = "";
#        $approxWindowRemaining = $approxWindow;
#      }
#    }
#  }
#
#  if(length($currentregion) > 0) {
#    push(@retVal, $currentregion);
#  }

  return \@retVal;
}

sub mergeFilterBCF {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $nbJobs = LoadConfig::getParam( $rH_cfg, 'mpileup', 'approxNbJobs' );
  my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);

  my $bcfDir = 'variants/rawBCF/';
  my $outputDir = 'variants/';

  my $command = SAMtools::mergeFilterBCF($rH_cfg, 'allSamples', $bcfDir, $outputDir, $rA_regions);
  my $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeFilterBCF", undef, 'MERGEBCF', $jobDependency, 'allSamples', $command);
  return '$'.$mergeJobId;
}

sub flagMappability {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $outputVCF = 'variants/allSamples.merged.flt.mil.vcf';
  my $inputVCF = 'variants/allSamples.merged.flt.vcf';
  my $command = VCFtools::annotateMappability($rH_cfg, $inputVCF, $outputVCF);
  my $milJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "flagMappability", undef, 'MAPPABILITY', $jobDependency, 'allSamples', $command);
  return '$'.$milJobId;
}

sub snpIDAnnotation {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $inputVCF = 'variants/allSamples.merged.flt.mil.vcf';
  my $vcfOutput = 'variants/allSamples.merged.flt.mil.snpId.vcf';

  my $command = SnpEff::annotateDbSnp($rH_cfg, $inputVCF, $vcfOutput);
  my $snpEffJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "snpIDAnnotation", undef, 'SNPID', $jobDependency, 'allSamples', $command);
 }


sub metricsSNV {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;
  
}


sub deliverable {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;


#   my $jobDependency ;
#   if($depends > 0) {
#     my $trimDependency = $globalDep{'trimMetrics'}{'trimMetrics'};
#     if (defined($trimDependency) && length($trimDependency) > 0) {
#       $jobDependency .= $trimDependency .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
#     }
#     my $alignDependency = $globalDep{'alignMetrics'}{'alignMetrics'};
#     if (defined($alignDependency) && length($alignDependency) > 0) {
#       $jobDependency .= $alignDependency .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
#     }
#     my $goDependency = $globalDep{'goseq'}{'goseq'};
#     if (defined($goDependency) && length($goDependency) > 0) {
#       $jobDependency .= $goDependency .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
#     }
#   }
# 
#   if (defined($jobDependency) && length($jobDependency) > 0) {
#     $jobDependency = substr $jobDependency, 0, -1 ;
#   }
# 
#   my $command = GqSeqUtils::clientReport($rH_cfg,  $configFile, $workDirectory) ;
# 
#   my $deliverableJobId = undef;
#   if(defined($command) && length($command) > 0) {
#     print "mkdir -p deliverable/output_jobs\n";
#     $deliverableJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "deliverable", 'REPORT', 'RNAREPORT', $jobDependency , undef, $command, 'deliverable' , $workDirectory);
#     $deliverableJobId = '$' .$deliverableJobId ;
#   }
# 
#   return $deliverableJobId;
}


1;
