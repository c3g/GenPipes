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

use BWA;
use GATK;
use IGVTools;
use LoadConfig;
use Picard;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SubmitToCluster;
use Trimmomatic;
#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'trimAndAlign'});
push(@steps, {'name' => 'mergeLanes'});
push(@steps, {'name' => 'indelRealigner'});
push(@steps, {'name' => 'mergeRealigned'});
push(@steps, {'name' => 'fixmate'});
push(@steps, {'name' => 'markDup'});
push(@steps, {'name' => 'metrics'});
push(@steps, {'name' => 'snpAndIndelBCF'});
push(@steps, {'name' => 'mergeFilterBCF'});
push(@steps, {'name' => 'filterNStretches'});
push(@steps, {'name' => 'flagMappability'});
push(@steps, {'name' => 'snpIDAnnotation'});
push(@steps, {'name' => 'snpEffect'});
push(@steps, {'name' => 'dbNSFPAnnotation'});
push(@steps, {'name' => 'indexVCF'});

#push(@steps, {'name' => 'crestSClip'});
push(@steps, {'name' => 'sortQname'});
#push(@steps, {'name' => 'countTelomere'});
push(@steps, {'name' => 'fullPileup'});
#push(@steps, {'name' => 'countTelomere'});
#  print "Step 14: filter N streches\n";
#  print "Step 15: flag mappability\n";
#  print "Step 16: snp annotation\n";
#  print "Step 17: snp effect prediction\n";
#  print "Step 18: gene descriptions and GO terms\n";
#  print "Strp 19: dbNSFP annotations\n";
#  print "Step 20: Cosmic annotations\n";


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
  for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};

    SubmitToCluster::initSubmit(\%cfg, $sampleName);
    for(my $current = $opts{'s'}-1; $current <= ($opts{'e'}-1); $current++) {
       my $fname = $steps[$current]->{'name'};
       my $subref = \&$fname;

       # Tests for the first step in the list. Used for dependencies.
       &$subref($current != ($opts{'s'}-1), \%cfg, $sampleName, $rAoH_sampleLanes, $rAoH_seqDictionary); 
    }
  }  
}

sub trimAndAlign {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  print "BWA_JOB_IDS=\"\"\n";
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo);
    my $trimJobIdVarName=undef;
    if(length($rH_trimDetails->{'command'}) > 0) {
      $trimJobIdVarName = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM', undef, $sampleName, $rH_trimDetails->{'command'});
      $trimJobIdVarName = '$'.$trimJobIdVarName;
      #TODO calcReadCounts.sh
    }

    my $rA_commands = BWA::aln($rH_cfg, $sampleName, $rH_laneInfo, $rH_trimDetails->{'pair1'}, $rH_trimDetails->{'pair2'}, $rH_trimDetails->{'single1'}, $rH_trimDetails->{'single2'});
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
  return '$BWA_JOB_IDS';
}

sub mergeLanes {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '$BWA_JOB_IDS';
  }

  my $command = Picard::merge($rH_cfg, $sampleName, $rAoH_sampleLanes);
  my $mergeJobId = undef;
  if(defined($command) && length($command) > 0) {
    $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "merge", undef, 'MERGELANES', $jobDependency, $sampleName, $command);
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
  my $outputPrefix = $sampleName.'/'.$sampleName.'.sorted.dup.all.coverage';
  $command = GATK::genomeCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  my $genomeCoverageJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "genomeCoverage", undef, 'GENOMECOVERAGE', $jobDependency, $sampleName, $command);

  # Compute CCDS coverage
  $outputPrefix = $sampleName.'/'.$sampleName.'.sorted.dup.CCDS.coverage';
  $command = GATK::targetCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  my $targetCoverageJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "targetCoverage", undef, 'TARGETCOVERAGE', $jobDependency, $sampleName, $command);

  # Generate IGV track
  $command = IGVTools::computeTDF($rH_cfg, $sampleName, $bamFile);
  my $igvtoolsTDFJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "computeTDF", undef, 'IGVTOOLS', $jobDependency, $sampleName, $command);

  # Compute flags
  my $output = $sampleName.'/'.$sampleName.'.sorted.dup.bam.flagstat';
  $command = SAMtools::flagstat($rH_cfg, $sampleName, $bamFile, $output);
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

sub snpAndIndelBCF {
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
  
  my $outputDir = LoadConfig::getParam($rH_cfg, "mpileup", 'sampleOutputRoot') . $sampleName."/rawBCF/";

  print 'mkdir -p '.$outputDir."\n";
  print "MPILEUP_JOB_IDS=\"\"\n";
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $snvWindow = LoadConfig::getParam($rH_cfg, 'mpileup', 'snvCallingWindow');
    my $seqName = $rH_seqInfo->{'name'};
    if($snvWindow ne "") {
      my $rA_regions = generateWindows($rH_seqInfo, $snvWindow);
      for my $region (@{$rA_regions}) {
        my $command = SAMtools::mpileup($rH_cfg, $sampleName, $bamFile, $region, $outputDir);
        my $mpileupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $region, 'MPILEUP', $jobDependency, $sampleName, $command);
        $mpileupJobId = '$'.$mpileupJobId;
        print 'MPILEUP_JOB_IDS=${MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$mpileupJobId."\n";
      }
    }
    else {
      my $command = SAMtools::mpileup($rH_cfg, $sampleName, $bamFile, $seqName, $outputDir);
      my $mpileupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $seqName, 'MPILEUP', $jobDependency, $sampleName, $command);
      $mpileupJobId = '$'.$mpileupJobId;
      print 'MPILEUP_JOB_IDS=${MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$mpileupJobId."\n";
    }
  }

  return '${MPILEUP_JOB_IDS}';
}

sub generateWindows {
  my $rH_seqInfo = shift;
  my $snvWindow = shift;

  my @retVal;
  for(my $idx=1; $idx <= $rH_seqInfo->{'size'}; $idx += $snvWindow) {
    my $end = $idx+$snvWindow-1;
    if($end > $rH_seqInfo->{'size'}) {
      $end = $rH_seqInfo->{'size'};
    }

    my $region = $rH_seqInfo->{'name'}.':'.$idx.'-'.$end;
    push(@retVal, $region);
  }

  return \@retVal;
}

sub mergeFilterBCF {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MPILEUP_JOB_IDS}';
  }

  my $snvWindow = LoadConfig::getParam($rH_cfg, 'mpileup', 'snvCallingWindow');

  my $bcfDir = LoadConfig::getParam($rH_cfg, "mergeFilterBCF", 'sampleOutputRoot') . $sampleName."/rawBCF/";
  my $outputDir = LoadConfig::getParam($rH_cfg, "mergeFilterBCF", 'sampleOutputRoot') . $sampleName.'/';

  my @seqNames;
  if($snvWindow ne "") {
    for my $rH_seqInfo (@$rAoH_seqDictionary) {
      my $rA_regions = generateWindows($rH_seqInfo, $snvWindow);
      push(@seqNames, @{$rA_regions});
    }
  }
  else {
    for my $rH_seqInfo (@$rAoH_seqDictionary) {
      push(@seqNames, $rH_seqInfo->{'name'});
    }
  }
  my $command = SAMtools::mergeFilterBCF($rH_cfg, $sampleName, $bcfDir, $outputDir, \@seqNames);
  my $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeFilterBCF", undef, 'MERGEBCF', $jobDependency, $sampleName, $command);
  return $mergeJobId;
}

1;
