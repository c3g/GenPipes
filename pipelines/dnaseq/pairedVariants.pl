#!/usr/bin/perl

=head1 NAME

I<pairedVariants>

=head1 SYNOPSIS

pairedVariants.pl

=head1 DESCRIPTION

B<pairedVariants> Pipeline to call snps, indels and Structural vairants.

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

#BEGIN{
#    #Makesure we can find the GetConfig::LoadModules module relative to this script install
#    use File::Basename;
#    use Cwd 'abs_path';
#    my ( undef, $mod_path, undef ) = fileparse( abs_path(__FILE__) );
#    unshift @INC, $mod_path."lib";
#
#}


# Dependencies
#--------------------
use Getopt::Std;

use LoadConfig;
use BAMtools;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SubmitToCluster;
use SVtools;
#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'snpAndIndelBCF'});
push(@steps, {'name' => 'mergeFilterBCF'});
push(@steps, {'name' => 'filterNStretches'});
push(@steps, {'name' => 'flagMappability'});
push(@steps, {'name' => 'snpAnnotation'});
push(@steps, {'name' => 'snpEffect'});
push(@steps, {'name' => 'geneDescriptions'});
push(@steps, {'name' => 'dbNSFP'});
push(@steps, {'name' => 'COSMIC'});
push(@steps, {'name' => 'DNAC'});
push(@steps, {'name' => 'filterDNAC'});
push(@steps, {'name' => 'Control-Freec'});

&main();

sub printUsage {
  print "\nUsage: perl ".$0." project.csv first_step last_step\n";
  print "\t-c  config file\n";
  print "\t-s  start step, inclusive\n";
  print "\t-e  end step, inclusive\n";
  print "\t-n  paired sample sheet. Format: 'SAMPLE,NORMAL,TUMOR'\n";
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
  my $rAoH_samplePairs = SampleSheet::parsePairedSampleSheet($opts{'n'});
  my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);

  my $latestBam;
  for my $rH_samplePair  (@$rAoH_samplePairs) {
    SubmitToCluster::initSubmit(\%cfg, $rH_samplePair->{'sample'});
    for(my $current = $opts{'s'}-1; $current <= ($opts{'e'}-1); $current++) {
       my $fname = $steps[$current]->{'name'};
       my $subref = \&$fname;

       # Tests for the first step in the list. Used for dependencies.
       &$subref($current != ($opts{'s'}-1), \%cfg, $rH_samplePair, $rAoH_seqDictionary); 
    }
  }  
}

sub snpAndIndelBCF {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.bam';
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.bam';
  my $outputDir = LoadConfig::getParam($rH_cfg, "mpileupPaired", 'sampleOutputRoot') . $sampleName."/rawBCF/";

  print 'mkdir -p '.$outputDir."\n";
  print "MPILEUP_JOB_IDS=\"\"\n";
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $snvWindow = LoadConfig::getParam($rH_cfg, 'mpileup', 'snvCallingWindow');
    my $seqName = $rH_seqInfo->{'name'};
    if($snvWindow ne "") {
      my $rA_regions = generateWindows($rH_seqInfo, $snvWindow);
      for my $region (@{$rA_regions}) {
        my $command = SAMtools::mpileupPaired($rH_cfg, $sampleName, $normalBam, $tumorBam, $region, $outputDir);
        my $mpileupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $region, 'MPILEUP', undef, $sampleName, $command);
        $mpileupJobId = '$'.$mpileupJobId;
        print 'MPILEUP_JOB_IDS=${MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$mpileupJobId."\n";
      } 
    }
    else {
      my $command = SAMtools::mpileupPaired($rH_cfg, $sampleName, $normalBam, $tumorBam, $seqName, $outputDir);
      my $mpileupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $seqName, 'MPILEUP', undef, $sampleName, $command);
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
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MPILEUP_JOB_IDS}';
  }

  my $snvWindow = LoadConfig::getParam($rH_cfg, 'mpileup', 'snvCallingWindow');

  my $sampleName = $rH_samplePair->{'sample'};
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

sub DNAC {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.bam';
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.bam';
  my $outputDir = LoadConfig::getParam($rH_cfg, "DNAC", 'sampleOutputRoot') . $sampleName.'/DNAC/';
  
  print 'mkdir -p '.$outputDir."\n";
  print "DNAC_JOB_IDS=\"\"\n";
  
  my $DNAC500File = $outputDir.$sampleName.'.DNAC_500';
  my $DNAC1000File = $outputDir.$sampleName.'.DNAC_1000.';
  my $DNAC30000File = $outputDir.$sampleName.'.DNAC_30000';
  my $bin500File = $DNAC500File.'.bins.tsv';
  my $bin1000File = $DNAC1000File.'.bins.tsv';
  my $bin30000File = $DNAC30000File.'.bins.tsv';

  my $command500 = BAMtools::countBins($rH_cfg, $sampleName, $tumorBam, 500, 'chr', $bin500File, $normalBam);
  my $command1000 = BAMtools::countBins($rH_cfg, $sampleName, $tumorBam, 1000, 'chr', $bin1000File, $normalBam);
  my $command30000 = BAMtools::countBins($rH_cfg, $sampleName, $tumorBam, 30000, 'genome', $bin30000File, $normalBam);

  $command500 .= ' && '.SVtools::runPairedDNAC($rH_cfg, $sampleName, $bin500File, $DNAC500File, 500);
  $command1000 .= ' && '.SVtools::runPairedDNAC($rH_cfg, $sampleName, $bin1000File, $DNAC1000File, 1000);
  $command30000 .= ' && '.SVtools::runPairedDNAC($rH_cfg, $sampleName, $bin30000File, $DNAC30000File, 30000);

  $command500 .= ' && '.SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC500File, $DNAC500File.'.bed', 1);
  $command1000 .= ' && '.SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC1000File, $DNAC1000File.'.bed', 1);
  $command30000 .= ' && '.SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC30000File, $DNAC30000File.'.bed', 2);

  my $dnac500JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC", undef, 'DNAC_500', $jobDependency, $sampleName, $command500);
  print 'DNAC_JOB_IDS=${DNAC_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$dnac500JobId."\n";
  my $dnac1000JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC", undef, 'DNAC_1000', $jobDependency, $sampleName, $command1000);
  print 'DNAC_JOB_IDS=${DNAC_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$dnac1000JobId."\n";
  my $dnac30000JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC", undef, 'DNAC_30000', $jobDependency, $sampleName, $command30000);
  print 'DNAC_JOB_IDS=${DNAC_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$dnac30000JobId."\n";

  return $dnac30000JobId;
}
1;
