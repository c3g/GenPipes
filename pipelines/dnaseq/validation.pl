#!/usr/bin/perl

=head1 NAME

I<validation>

=head1 SYNOPSIS

validation.pl

=head1 DESCRIPTION

B<validation> Is the SNV/Indel validation pipeline

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
use LoadConfig;
use Picard;
use SampleSheet;
use SequenceDictionaryParser;
use SubmitToCluster;
use Trimmomatic;
#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'trim'});
push(@steps, {'name' => 'align'});
push(@steps, {'name' => 'metrics'});

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
  #my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);
  my $rAoH_seqDictionary = undef;

  SubmitToCluster::initPipeline;

  my $latestBam;
  for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};

    for my $rH_laneInfo (@$rAoH_sampleLanes) {
      for(my $current = $opts{'s'}-1; $current <= ($opts{'e'}-1); $current++) {
        my $fname = $steps[$current]->{'name'};
        my $subref = \&$fname;
        # Tests for the first step in the list. Used for dependencies.
        &$subref($current != ($opts{'s'}-1), \%cfg, $sampleName, $rH_laneInfo, $rAoH_seqDictionary);
      } 
    }
  }  
}

sub trim {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rH_laneInfo  = shift;
  my $rAoH_seqDictionary = shift;

  my $outputDir = 'reads/'.$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
  print 'mkdir -p '.$outputDir."\n";

  my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo, $outputDir);
  my $trimJobIdVarName=undef;
  if(length($rH_trimDetails->{'command'}) > 0) {
    $trimJobIdVarName = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM', undef, $sampleName, $rH_trimDetails->{'command'});
  }
  return $trimJobIdVarName;
}

sub align {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rH_laneInfo  = shift;
  my $rAoH_seqDictionary = shift;

  my $minQuality = $rH_cfg->{'trim.minQuality'};
  my $minLength = $rH_cfg->{'trim.minLength'};
  my $adapterFile = $rH_cfg->{'trim.adapterFile'};

  my $inputDir = 'reads/'.$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
  my $outputDir = 'alignment/'.$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
  my $inputFastqPair1Name = $inputDir . $sampleName.'.t'.$minQuality.'l'.$minLength.'.pair1.fastq.gz';
  my $inputFastqPair2Name = $inputDir . $sampleName.'.t'.$minQuality.'l'.$minLength.'.pair2.fastq.gz';
  my $inputFastqSingleName = $inputDir . $sampleName.'.t'.$minQuality.'l'.$minLength.'.single1.fastq.gz';
  my $outputPrefix = $outputDir . $sampleName;
  print 'mkdir -p '.$outputDir."\n";

  my $jobDep = "";
  print "BWA_JOB_IDS=\"\"\n";
  my $trimJobIdVarName = undef;
  if($depends) {
    $trimJobIdVarName = '$TRIM_JOB_ID';
  }

  my $rgId = $rH_laneInfo->{'libraryBarcode'} . "_" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
  my $rgSampleName = $rH_laneInfo->{'name'};
  my $rgLibrary = $rH_laneInfo->{'libraryBarcode'};
  my $rgPlatformUnit = $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
  my $rgCenter = LoadConfig::getParam( $rH_cfg, 'mem', 'bwaInstitution' );

  if($rH_laneInfo->{'runType'} eq "PAIRED_END") {
    my $rA_commands = BWA::mem($rH_cfg, $sampleName, $inputFastqPair1Name, $inputFastqPair2Name, undef, $outputPrefix.'.paired', $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter);
    my $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mem", 'memPair.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA', $trimJobIdVarName, $sampleName, $rA_commands->[0]);
    $bwaJobId = '$'.$bwaJobId;
    print 'BWA_JOB_IDS='.$bwaJobId.LoadConfig::getParam($rH_cfg, 'mem', 'clusterDependencySep').'${BWA_JOB_IDS}'."\n";
    
    # fake it and take the single1 end
    $rA_commands = BWA::mem($rH_cfg, $sampleName, undef, undef, $inputFastqSingleName, $outputPrefix.'.single', $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter);
    $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mem", 'memSingle.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA', $trimJobIdVarName, $sampleName, $rA_commands->[0]);
    $bwaJobId = '$'.$bwaJobId;
    print 'BWA_JOB_IDS='.$bwaJobId.LoadConfig::getParam($rH_cfg, 'mem', 'clusterDependencySep').'${BWA_JOB_IDS}'."\n";

    my $outputPairedBAM = $outputDir . $sampleName.'.paired.sorted.bam';
    my $outputSingleBAM = $outputDir . $sampleName.'.single.sorted.bam';
    my $outputBAM = $outputDir . $sampleName.'.sorted.bam';
    my @inputBams = ($outputPairedBAM, $outputSingleBAM);
    my $command = Picard::mergeFiles($rH_cfg, $sampleName, \@inputBams, $outputBAM);
    my $mergeJobId = undef;
    if(defined($command) && length($command) > 0) {
      $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergePairs", undef, 'MERGEPAIRS', '$BWA_JOB_IDS', $sampleName, $command);
    }
    if(defined($mergeJobId)) {
      $jobDep = '$'.$mergeJobId;
    }
  }
  else {
    my $rA_commands = BWA::mem($rH_cfg, $sampleName, $inputFastqPair1Name, $inputFastqPair2Name, $inputFastqSingleName, $outputPrefix, $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter);
    my $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mem", 'memSingle.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA', $trimJobIdVarName, $sampleName, $rA_commands->[0]);
    $bwaJobId = '$'.$bwaJobId;
    print 'BWA_JOB_IDS='.$bwaJobId.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').'${BWA_JOB_IDS}'."\n";
    $jobDep = '$BWA_JOB_IDS';
  }

  return $jobDep;
}

sub metrics {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rH_laneInfo  = shift;
  my $rAoH_seqDictionary = shift;

  my $bwaJobId = undef;
  if($depends) {
    $bwaJobId = '$MERGEPAIRS_JOB_ID';
  }
  my $laneDirectory = 'alignment/'.$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
  # Compute target coverage
  my $command = GATK::targetCoverage($rH_cfg, $sampleName, $laneDirectory . $sampleName.'.sorted.bam', $laneDirectory . $sampleName.'.sorted.targetCoverage');
  my $genomeCoverageJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "targetCoverage", undef, 'TARGETCOVERAGE', $bwaJobId, $sampleName, $command);
}

1;
