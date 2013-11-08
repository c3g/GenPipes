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

BEGIN {
  # Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
  use File::Basename;
  use Cwd 'abs_path';
  my (undef, $mod_path, undef) = fileparse(abs_path(__FILE__));
  unshift @INC, $mod_path . "../../lib";
}


# Dependencies
#--------------------
use Getopt::Std;
use POSIX;

use LoadConfig;
use BAMtools;
use Breakdancer;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SnpEff;
use SubmitToCluster;
use SVtools;
use Tools;
use Version;
use VCFtools;
use Pindel;
use Cfreec;

#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'snpAndIndelBCF', 'stepLoop' => 'sample', 'parentStep' => undef});
push(@steps, {'name' => 'mergeFilterBCF', 'stepLoop' => 'sample', 'parentStep' => 'snpAndIndelBCF'});
push(@steps, {'name' => 'filterNStretches', 'stepLoop' => 'sample', 'parentStep' => 'mergeFilterBCF'});
push(@steps, {'name' => 'flagMappability', 'stepLoop' => 'sample', 'parentStep' => 'filterNStretches'});
push(@steps, {'name' => 'snpIDAnnotation', 'stepLoop' => 'sample', 'parentStep' => 'flagMappability'});
push(@steps, {'name' => 'snpEffect', 'stepLoop' => 'sample', 'parentStep' => 'snpIDAnnotation'});
push(@steps, {'name' => 'dbNSFPAnnotation', 'stepLoop' => 'sample', 'parentStep' => 'snpEffect'});
push(@steps, {'name' => 'DNAC', 'stepLoop' => 'sample', 'parentStep' => undef});
push(@steps, {'name' => 'Breakdancer', 'stepLoop' => 'sample', 'parentStep' => undef});
push(@steps, {'name' => 'Pindel', 'stepLoop' => 'sample', 'parentStep' => undef});
push(@steps, {'name' => 'ControlFreec', 'stepLoop' => 'sample', 'parentStep' => undef});

my %globalDep;
for my $stepName (@steps) {
  $globalDep{$stepName -> {'name'} } ={};
}

&main();

sub printUsage {
  print "Version: ".$Version::version."\n";
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

  SubmitToCluster::initPipeline;

  my $latestBam;
  for my $rH_samplePair  (@$rAoH_samplePairs) {
    for(my $currentStep = $opts{'s'}-1; $currentStep <= ($opts{'e'}-1); $currentStep++) {
      my $fname = $steps[$currentStep]->{'name'};
      my $subref = \&$fname;

      # Tests for the first step in the list. Used for dependencies.
      my $jobIdVar = &$subref($currentStep, \%cfg, $rH_samplePair, $rAoH_seqDictionary); 
      $globalDep{$fname}->{$rH_samplePair->{'sample'}} = $jobIdVar;
    }
  }  
}

sub snpAndIndelBCF {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = 'alignment/'.$rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.recal.bam';
  my $tumorBam = 'alignment/'.$rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.recal.bam';
  my $outputDir = 'pairedVariants/' . $sampleName."/rawBCF/";

  print 'mkdir -p '.$outputDir."\n";
  print "MPILEUP_JOB_IDS=\"\"\n";

  my $nbJobs = LoadConfig::getParam( $rH_cfg, 'mpileup', 'approxNbJobs' );
  my $jobId;
  if (defined($nbJobs) && $nbJobs > 1) {
    my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);
    for my $region (@{$rA_regions}) {
      my $rO_job = SAMtools::mpileupPaired($rH_cfg, $sampleName, $normalBam, $tumorBam, $region, $outputDir);
      if(!$rO_job->isUp2Date()) {
        $region =~ s/:/_/g;
        SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $region, 'MPILEUP', $jobDependency, $sampleName, $rO_job);
        if(!defined($jobId)) {
          $jobId = '${MPILEUP_JOB_IDS}';
          print 'MPILEUP_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
        }
        else {
          print 'MPILEUP_JOB_IDS=${MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
        }
      }
    }
  }
  else {
    my $rO_job = SAMtools::mpileupPaired($rH_cfg, $sampleName, $normalBam, $tumorBam, undef, $outputDir);
    SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", undef, 'MPILEUP', undef, $sampleName, $rO_job);
    if(!defined($jobId)) {
      $jobId = '${MPILEUP_JOB_IDS}';
      print 'MPILEUP_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
    }
    else {
      print 'MPILEUP_JOB_IDS=${MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
    }
  }

  return $jobId;
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
  return \@retVal;
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
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $nbJobs = LoadConfig::getParam( $rH_cfg, 'mpileup', 'approxNbJobs' );
  my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);

  my $sampleName = $rH_samplePair->{'sample'};
  my $bcfDir = 'pairedVariants/' . $sampleName."/rawBCF/";
  my $outputDir = 'pairedVariants/' . $sampleName.'/'; 

  my $rO_job = SAMtools::mergeFilterBCF($rH_cfg, $sampleName, $bcfDir, $outputDir, $rA_regions);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeFilterBCF", undef, 'MERGEBCF', $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub filterNStretches {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $inputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.vcf';
  my $outputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.NFiltered.vcf';

  my $rO_job = Tools::filterNStretches($rH_cfg, $sampleName, $inputVCF, $outputVCF);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "filterNStretches", undef, 'FILTERN', $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub flagMappability {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $inputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.NFiltered.vcf';
  my $outputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.mil.vcf';

  # Use mergeFilterBCF to make sure we have the right path
  my $rO_job = VCFtools::annotateMappability($rH_cfg, $inputVCF, $outputVCF);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "flagMappability", undef, 'MAPPABILITY', $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub snpIDAnnotation {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $inputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.mil.vcf';
  my $vcfOutput = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.mil.snpId.vcf';
  my $rO_job = SnpEff::annotateDbSnp($rH_cfg, $inputVCF, $vcfOutput);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "snpIDAnnotation", undef, 'SNPID', $jobDependency, $sampleName, $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

sub snpEffect {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $inputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.mil.snpId.vcf';
  my $vcfOutput = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.mil.snpId.snpeff.vcf';


  my $rO_job = SnpEff::computeEffects($rH_cfg, $inputVCF, $vcfOutput, 1);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "snpEffect", undef, 'SNPEFF', $jobDependency, $sampleName, $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

sub dbNSFPAnnotation {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $inputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.mil.snpId.snpeff.vcf';
  my $vcfOutput = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.mil.snpId.snpeff.dbnsfp.vcf';
  my $rO_job = SnpEff::annotateDbNSFP($rH_cfg, $inputVCF, $vcfOutput);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "dbNSFPAnnotation", undef, 'DBNSFP', $jobDependency, 'allSamples', $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

sub DNAC {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.bam';
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.bam';
  my $outputDir = 'pairedVariants/' . $sampleName.'/DNAC/';
  
  print 'mkdir -p '.$outputDir."\n";
  print "DNAC_JOB_IDS=\"\"\n";
  
  my $DNAC500File = $outputDir.$sampleName.'.DNAC_500';
  my $DNAC1000File = $outputDir.$sampleName.'.DNAC_1000';
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

  $command500 .= ' && '.SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC500File.'.txt', $DNAC500File.'.filteredSV', 1);
  $command1000 .= ' && '.SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC1000File.'.txt', $DNAC1000File.'.filteredSV', 1);
  $command30000 .= ' && '.SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC30000File.'.txt', $DNAC30000File.'.filteredSV', 2);

  my $dnac500JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC_500", undef, 'DNAC_500', $jobDependency, $sampleName, $command500);
  print 'DNAC_JOB_IDS=${DNAC_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$dnac500JobId."\n";
  my $dnac1000JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC_1000", undef, 'DNAC_1000', $jobDependency, $sampleName, $command1000);
  print 'DNAC_JOB_IDS=${DNAC_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$dnac1000JobId."\n";
  my $dnac30000JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC_30000", undef, 'DNAC_30000', $jobDependency, $sampleName, $command30000);
  print 'DNAC_JOB_IDS=${DNAC_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$dnac30000JobId."\n";

  return $dnac30000JobId;
}

sub Breakdancer {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.bam';
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.bam';
  my $outputDir = 'pairedVariants/' . $sampleName.'/breakdancer/';

  print 'mkdir -p '.$outputDir."\n";
  print "BRD_JOB_IDS=\"\"\n";

  my $normalOutput = $outputDir.$sampleName.'.brdN.cfg';
  my $tumorOutput = $outputDir.$sampleName.'.brdT.cfg';
  my $sampleCFGOutput = $outputDir.$sampleName.'.brd.cfg';
  my $command = Breakdancer::bam2cfg($rH_cfg,$sampleName,$normalBam,$normalOutput, LoadConfig::getParam($rH_cfg, 'Breakdancer', 'normalStdDevCutoff'));
  $command .= ' & '.Breakdancer::bam2cfg($rH_cfg,$sampleName,$tumorBam,$tumorOutput, LoadConfig::getParam($rH_cfg, 'Breakdancer', 'tumorStdDevCutoff'));
  $command .= ' & wait && cat '.$normalOutput.' '.$tumorOutput. ' > '.$sampleCFGOutput;
  my $brdCFGJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "Breakdancer", undef, 'BRD_CFG', $jobDependency, $sampleName, $command);
  $brdCFGJobId = '$'.$brdCFGJobId;

  my $outputTRPrefix = $outputDir.$sampleName.'.brd.TR';
  $command = Breakdancer::pairedBRDITX($rH_cfg,$sampleName,$sampleCFGOutput,$outputTRPrefix);
  my $brdTRJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "Breakdancer", 'TR', 'BRD_TR', $brdCFGJobId, $sampleName, $command);
  $brdTRJobId = '$'.$brdTRJobId;

  print "BRD_JOB_IDS=\"".$brdTRJobId."\"\n";
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};

    my $outputPrefix = $outputDir.$sampleName.'.brd.'.$seqName;
    $command = Breakdancer::pairedBRD($rH_cfg,$sampleName,$seqName,$sampleCFGOutput,$outputPrefix);
    my $brdJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "Breakdancer", $seqName, 'BRD', $brdCFGJobId, $sampleName, $command);
    $brdJobId = '$'.$brdJobId;

    print 'BRD_JOB_IDS=${BRD_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$brdJobId."\n";
  }
  ## merge results
  my $outputPrefix = $outputDir.$sampleName.'.brd';
  $command = Breakdancer::mergeCTX($rH_cfg,$outputPrefix) ;
  ##filter results
  my $brdCallsFile = $outputPrefix .'.ctx';
  $command .= ' && '.SVtools::filterBrD($rH_cfg, $sampleName, $brdCallsFile, $outputPrefix.'.filteredSV', $normalBam, $tumorBam);
  my $brdMergeFilterJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "filterSV", $sampleName, 'FILTBRD', '${BRD_JOB_IDS}', $sampleName, $command);
  print 'FILTBRD_JOB_IDS=$'.$brdMergeFilterJobId."\n";
  return '$FILTBRD_JOB_IDS';
}

sub Pindel {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.bam';
  my $normalMetrics = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.all.metrics.insert_size_metrics';
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.bam';
  my $tumorMetrics = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.all.metrics.insert_size_metrics';
  my $outputDir = 'pairedVariants/' . $sampleName.'/pindel/';

  print 'mkdir -p '.$outputDir."\n";
  print "PI_JOB_IDS=\"\"\n";

  my $sampleCFGOutput = $outputDir.$sampleName.'.Pindel.conf';
  my $command = Pindel::pairedConfigFile($rH_cfg, $tumorMetrics, $normalMetrics, $tumorBam, $normalBam, $sampleCFGOutput);
  my $piCFGJobId; 
  if(defined($command) && length($command) > 0) {
    $piCFGJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "PindelCfg", undef, 'PI_CFG', $jobDependency, $sampleName, $command);
    $piCFGJobId = '$'.$piCFGJobId;
  }
  my $pi2filterDep;
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $chrFile = LoadConfig::getParam($rH_cfg, "Pindel", 'referenceGenomeByChromosome') .'/chr' .$seqName.'.fa';
    my $Brdresfile = 'pairedVariants/' . $sampleName.'/breakdancer/'.$sampleName.'.brd.'.$seqName.'ctx';
    my $outputPrefix = $outputDir .$sampleName .'.' .$seqName;
    my $outputTest= $outputDir .$sampleName ;
    my $BrdOption = '' ;
    if (-e $Brdresfile) {
      my $outputBrdPi = $outputDir .$sampleName .'.PI_BrD.calls.txt';
      $BrdOption .= ' -Q ' .$outputBrdPi .' -b ' .$Brdresfile  ;
    } 
    $command = Pindel::pairedPI($rH_cfg, $chrFile, $sampleCFGOutput, $outputPrefix, $outputTest, $BrdOption);
    my $piJobId ; 
    if(defined($command) && length($command) > 0) {
      $piJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "pindel", $seqName, 'PI', $piCFGJobId, $sampleName, $command);
      $piJobId = '$'.$piJobId;
      print 'PI_JOB_IDS=${PI_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$piJobId."\n";
      $pi2filterDep='${PI_JOB_IDS}';
    }
  }
  ## merge results
  my $outputPrefix = $outputDir.$sampleName;
  $command = Pindel::mergeChro($rH_cfg,$outputPrefix) ;
  ##filter results
  if(defined($command) && length($command) > 0) {
     $command .= ' && ';
  }
  $command .= SVtools::filterPI($rH_cfg, $sampleName, $outputPrefix, $outputPrefix.'.filteredSV', $normalBam, $tumorBam);
  my $piMergeFilterJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "filterSV", $sampleName .'PI', 'FILTPI', $pi2filterDep, $sampleName, $command);
  print 'FILTPI_JOB_IDS=$'.$piMergeFilterJobId."\n";
  return '$FILTPI_JOB_IDS';
}

sub ControlFreec {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.'.LoadConfig::getParam($rH_cfg, "ControlFreec", 'inputExtension');
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.'.LoadConfig::getParam($rH_cfg, "ControlFreec", 'inputExtension');
  my $outputDir = 'pairedVariants/' . $sampleName.'/controlFREEC/';
  my $sampleConfigFile = $outputDir.'/'.$sampleName.'.freec.cfg';
  

  print 'mkdir -p '.$outputDir."\n";
  print "CONTROL_FREEC_JOB_IDS=\"\"\n";

  my $command = Cfreec::pairedFreec($rH_cfg, $tumorBam, $normalBam, $sampleConfigFile, $outputDir);
  my $cFreecJobId; 
  if(defined($command) && length($command) > 0) {
    $cFreecJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "ControlFreec", undef, 'CONTROL_FREEC', $jobDependency, $sampleName, $command);
    $cFreecJobId = '$'.$cFreecJobId;
  }
  return $cFreecJobId;

}

1;

