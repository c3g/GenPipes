#!/usr/bin/env perl

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

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin/../../lib";

# Dependencies
#--------------------
use Getopt::Std;
use POSIX;

use LoadConfig;
use BVATools;
use Breakdancer;
use GATK;
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
push(@steps, {'name' => 'mutect', 'stepLoop' => 'sample', 'parentStep' => undef});
push(@steps, {'name' => 'mergeMuTect', 'stepLoop' => 'sample', 'parentStep' => 'mutect'});

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

  my $nbJobs = LoadConfig::getParam($rH_cfg, 'mpileup', 'approxNbJobs', 0, 'int');
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
    SubmitToCluster::printSubmitCmd($rH_cfg, "dbNSFPAnnotation", undef, 'DBNSFP', $jobDependency, $sampleName, $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

sub DNAC {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

#   my $jobDependency = undef;
#   my $parentStep = $steps[$stepId]->{'parentStep'};
#   if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
#     $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
#   }
# 
#   my $sampleName = $rH_samplePair->{'sample'};
#   my $inputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.NFiltered.vcf';
#   my $outputVCF = 'pairedVariants/' . $sampleName.'/'.$sampleName.'.merged.flt.mil.vcf';
# 
#   # Use mergeFilterBCF to make sure we have the right path
#   my $rO_job = VCFtools::annotateMappability($rH_cfg, $inputVCF, $outputVCF);
#   if(!$rO_job->isUp2Date()) {
#     SubmitToCluster::printSubmitCmd($rH_cfg, "flagMappability", undef, 'MAPPABILITY', $jobDependency, $sampleName, $rO_job);
#   }
#   return $rO_job->getCommandJobId(0);
#   my $trimDependency = $ro_trimJob->getCommandJobId(0);
   my $sampleName = $rH_samplePair->{'sample'};
  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($parentStep) && defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $normalBam = 'alignment/'.$rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.recal.bam';
  my $tumorBam = 'alignment/'.$rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.recal.bam';
  my $outputDir = 'pairedVariants/' . $sampleName.'/DNAC/';
  
  print 'mkdir -p '.$outputDir."\n";
  print "DNAC_JOB_IDS=\"\"\n";
  
  my $DNAC500File = $outputDir.$sampleName.'.DNAC_500';
  my $DNAC1000File = $outputDir.$sampleName.'.DNAC_1000';
  my $DNAC30000File = $outputDir.$sampleName.'.DNAC_30000';
  my $bin500File = $DNAC500File.'.bins.tsv';
  my $bin1000File = $DNAC1000File.'.bins.tsv';
  my $bin30000File = $DNAC30000File.'.bins.tsv';
  ## run binnng step with BVA
  my $rO_binCout500_job = BVATools::countBins($rH_cfg, $sampleName, $tumorBam, 500, 'chr', $bin500File, $normalBam);
  if(!$rO_binCout500_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC_BIN", 'BIN_COUNT_500', 'BIN_COUNT_500', $jobDependency, $sampleName, $rO_binCout500_job);
  }
  my $rO_binCout1000_job = BVATools::countBins($rH_cfg, $sampleName, $tumorBam, 1000, 'chr', $bin1000File, $normalBam);
  if(!$rO_binCout1000_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC_BIN", 'BIN_COUNT_1000', 'BIN_COUNT_1000', $jobDependency, $sampleName, $rO_binCout1000_job);
  }
  my $rO_binCout30000_job = BVATools::countBins($rH_cfg, $sampleName, $tumorBam, 30000, 'genome', $bin30000File, $normalBam);
  if(!$rO_binCout30000_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC_BIN", 'BIN_COUNT_30000', 'BIN_COUNT_30000', $jobDependency, $sampleName, $rO_binCout30000_job);
  }
   ## run CNV calling step with DNACRD
  my $GCmapFile500 =  LoadConfig::getParam($rH_cfg, 'DNAC', 'GCmapFileBasename') .'_bin500bp_GCMAP.bed';
  my $rO_dnac500_job = SVtools::runPairedDNAC($rH_cfg, $sampleName, $bin500File, $DNAC500File, 500, $GCmapFile500);
  if(!$rO_dnac500_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC", 'DNAC_500', 'DNAC_500', $rO_binCout500_job->getCommandJobId(0), $sampleName, $rO_dnac500_job);
  }
  my $GCmapFile1000 =  LoadConfig::getParam($rH_cfg, 'DNAC', 'GCmapFileBasename') .'_bin1kb_GCMAP.bed';
  my $rO_dnac1000_job = SVtools::runPairedDNAC($rH_cfg, $sampleName, $bin1000File, $DNAC1000File, 1000, $GCmapFile1000);
  if(!$rO_dnac1000_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC", 'DNAC_1000', 'DNAC_1000', $rO_binCout1000_job->getCommandJobId(0), $sampleName, $rO_dnac1000_job);
  }
  my $GCmapFile30000 =  LoadConfig::getParam($rH_cfg, 'DNAC', 'GCmapFileBasename') .'_bin30kb_GCMAP.bed';
  my $rO_dnac30000_job = SVtools::runPairedDNAC($rH_cfg, $sampleName, $bin30000File, $DNAC30000File, 30000,$GCmapFile30000);
  if(!$rO_dnac30000_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC", 'DNAC_30000', 'DNAC_30000', $rO_binCout30000_job->getCommandJobId(0), $sampleName, $rO_dnac30000_job);
  }
  ## run 1st filtering step with SVtools script
  my $rO_filter500_job = SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC500File.'_CNVcalls.txt', $DNAC500File.'.filteredSV', 1);
  if(!$rO_filter500_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC_FILT", 'FILTER1_DNAC_500', 'FILTER1_DNAC_500', $rO_dnac500_job->getCommandJobId(0), $sampleName, $rO_filter500_job);
  }
  my $rO_filter1000_job = SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC1000File.'_CNVcalls.txt', $DNAC1000File.'.filteredSV', 1);
  if(!$rO_filter1000_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "DNAC_FILT", 'FILTER1_DNAC_1000', 'FILTER1_DNAC_1000', $rO_dnac1000_job->getCommandJobId(0), $sampleName, $rO_filter1000_job);
  }
  my $rO_filter30000_job = SVtools::filterDNAC($rH_cfg, $sampleName, $DNAC30000File.'_CNVcalls.txt', $DNAC30000File.'.filteredSV', 2);
  if(!$rO_filter30000_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "filterSV", 'FILTER1_DNAC_30000', 'FILTER1_DNAC_30000', $rO_dnac30000_job->getCommandJobId(0), $sampleName, $rO_filter30000_job);
  }
  
  return $rO_filter500_job->getCommandJobId(0) .':' .$rO_filter1000_job->getCommandJobId(0) .':' .$rO_filter30000_job->getCommandJobId(0) ;
}

sub Breakdancer {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $sampleName = $rH_samplePair->{'sample'};
  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($parentStep) && defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

#  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = 'alignment/'.$rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.recal.bam';
  my $tumorBam = 'alignment/'.$rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.recal.bam';
  my $outputDir = 'pairedVariants/' . $sampleName.'/breakdancer/';

  print 'mkdir -p '.$outputDir."\n";
  print "BRD_JOB_IDS=\"\"\n";

  ##generates config files (tumor and normal)

  my $normalOutput = $outputDir.$sampleName.'.brdN.cfg';
  my $tumorOutput = $outputDir.$sampleName.'.brdT.cfg';
  my $sampleCFGOutput = $outputDir.$sampleName.'.brd.cfg';
  my $rO_normal_job = Breakdancer::bam2cfg($rH_cfg,$sampleName,$normalBam,$normalOutput, LoadConfig::getParam($rH_cfg, 'BRD', 'normalStdDevCutoff'));
  if(!$rO_normal_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "BRD_CONF",'CONF_N', 'CONF_N', $jobDependency, $sampleName, $rO_normal_job);
  }
  my $rO_tumor_job = Breakdancer::bam2cfg($rH_cfg,$sampleName,$tumorBam,$tumorOutput, LoadConfig::getParam($rH_cfg, 'BRD', 'tumorStdDevCutoff'));
  if(!$rO_tumor_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "BRD_CONF", 'CONF_T', 'CONF_T', $jobDependency, $sampleName, $rO_tumor_job);
  }

  ## merge tumor and normal conf output
  my $rO_mergeConf_job = Breakdancer::mergeConf($rH_cfg, $normalOutput, $tumorOutput, $sampleCFGOutput);
  if(!$rO_mergeConf_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "BRD_CONF", 'CONF_MERGE', 'CONF_MERGE',  $rO_normal_job->getCommandJobId(0) .LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_tumor_job->getCommandJobId(0) , $sampleName, $rO_mergeConf_job);
  }
  
  ##run BRD translocation
  my $outputTRPrefix = $outputDir.$sampleName.'.brd.TR';
  my $rO_brdT_job = Breakdancer::pairedBRDITX($rH_cfg,$sampleName,$sampleCFGOutput,$outputTRPrefix);
  if(!$rO_brdT_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "BRD", 'BRD_TR', 'BRD_TR',  $rO_mergeConf_job->getCommandJobId(0), $sampleName, $rO_brdT_job);
  }

  my $brdjobId = $rO_brdT_job->getCommandJobId(0);

  ##run BRD by chromosome
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $outputPrefix = $outputDir.$sampleName.'.brd.'.$seqName;
    my $rO_brdChro_job = Breakdancer::pairedBRD($rH_cfg,$sampleName,$seqName,$sampleCFGOutput,$outputPrefix);
    if(!$rO_brdChro_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "BRD", 'BRD_'.$seqName, 'BRD_'.$seqName,  $rO_mergeConf_job->getCommandJobId(0), $sampleName, $rO_brdChro_job);
      $brdjobId .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_brdChro_job->getCommandJobId(0);
    }
  }
  ## merge results
  my $outputPrefix = $outputDir.$sampleName.'.brd';
  my $rO_brdMerge_job = Breakdancer::mergeCTX($rH_cfg,$outputPrefix) ;
  if(!$rO_brdMerge_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "BRD_MERGE", 'BRD_MERGE', 'BRD_MERGE',  $brdjobId, $sampleName, $rO_brdMerge_job);
  }
  ##filter results
  my $brdCallsFile = $outputPrefix .'.ctx';
  my $rO_brdFilter_job = SVtools::filterBrD($rH_cfg, $sampleName, $brdCallsFile, $outputPrefix.'.filteredSV', $normalBam, $tumorBam);
  if(!$rO_brdFilter_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "filterSV", 'BRD_FILTER', 'BRD_FILTER',  $rO_brdMerge_job->getCommandJobId(0), $sampleName, $rO_brdFilter_job);
  }
  return $rO_brdFilter_job->getCommandJobId(0);
}

sub Pindel {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $sampleName = $rH_samplePair->{'sample'};
  my $jobDependency = undef;
#  my $parentStep = $steps[$stepId]->{'parentStep'};
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($parentStep) && defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

#  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = 'alignment/'.$rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.recal.bam';
  my $normalMetrics = 'alignment/'.$rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.isize.txt';
  my $tumorBam = 'alignment/'.$rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.recal.bam';
  my $tumorMetrics = 'alignment/'.$rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.isize.txt';
  my $outputDir = 'pairedVariants/' . $sampleName.'/pindel/';

  print 'mkdir -p '.$outputDir."\n";
  my $pi2filterDep ;
  ### run config
  my $sampleCFGOutput = $outputDir.$sampleName.'.Pindel.conf';
  my $rO_pindelConf_job = Pindel::pairedConfigFile($rH_cfg, $tumorMetrics, $normalMetrics, $tumorBam, $normalBam, $sampleCFGOutput);
  if(!$rO_pindelConf_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "PindelCfg", 'PI_CFG', 'PI_CFG',  $jobDependency, $sampleName, $rO_pindelConf_job);
      $pi2filterDep = $rO_pindelConf_job->getCommandJobId(0);
  }
  
   ### run pindel
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
    my $rO_pindelChro_job = Pindel::pairedPI($rH_cfg, $chrFile, $sampleCFGOutput, $outputPrefix, $outputTest, $BrdOption);
    if(!$rO_pindelChro_job->isUp2Date()) {
        SubmitToCluster::printSubmitCmd($rH_cfg, "Pindel", 'PI_'.$seqName, 'PI_'.$seqName,  $rO_pindelConf_job->getCommandJobId(0), $sampleName, $rO_pindelChro_job);
        $pi2filterDep .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_pindelChro_job->getCommandJobId(0);
    }
  }
  ## merge results
  my $outputPrefix = $outputDir.$sampleName;
  my $rO_pindelMerge_job = Pindel::mergeChro($rH_cfg,$outputPrefix) ;
  if(!$rO_pindelMerge_job->isUp2Date()) {
        SubmitToCluster::printSubmitCmd($rH_cfg, "Pindel_merge", 'PI_MERGE', 'PI_MERGE',  $pi2filterDep, $sampleName, $rO_pindelMerge_job);
  }
  ##filter results
  my $rO_pindelFilter_job= SVtools::filterPI($rH_cfg, $sampleName, $outputPrefix, $outputPrefix.'.filteredSV', $normalBam, $tumorBam);
  if(!$rO_pindelFilter_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "filterSV", 'FILTPI', 'FILTPI',$rO_pindelMerge_job->getCommandJobId(0) , $sampleName,$rO_pindelFilter_job);
  }
  return $rO_pindelFilter_job->getCommandJobId(0);
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

sub mutect {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = 'alignment/'.$rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.recal.bam';
  my $tumorBam = 'alignment/'.$rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.recal.bam';
  my $outputDir = 'pairedVariants/' . $sampleName."/rawMuTect/";

  print 'mkdir -p '.$outputDir."\n";
  print "MUTECT_JOB_IDS=\"\"\n";

  my $nbJobs = LoadConfig::getParam($rH_cfg, 'mutect', 'approxNbJobs', 0, 'int');
  my $jobId;
  if (defined($nbJobs) && $nbJobs > 1) {
    my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);
    for my $region (@{$rA_regions}) {
      my $rO_job = GATK::mutect($rH_cfg, $sampleName, $normalBam, $tumorBam, $region, $outputDir);
      if(!$rO_job->isUp2Date()) {
        $region =~ s/:/_/g;
        SubmitToCluster::printSubmitCmd($rH_cfg, "mutect", $region, 'MUTECT', $jobDependency, $sampleName, $rO_job);
        if(!defined($jobId)) {
          $jobId = '${MUTECT_JOB_IDS}';
          print 'MUTECT_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
        }
        else {
          print 'MUTECT_JOB_IDS=${MUTECT_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
        }
      }
    }
  }
  else {
    my $rO_job = GATK::mutect($rH_cfg, $sampleName, $normalBam, $tumorBam, undef, $outputDir);
    SubmitToCluster::printSubmitCmd($rH_cfg, "mutect", undef, 'MUTECT', undef, $sampleName, $rO_job);
    if(!defined($jobId)) {
      $jobId = '${MUTECT_JOB_IDS}';
      print 'MUTECT_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
    }
    else {
      print 'MUTECT_JOB_IDS=${MUTECT_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
    }
  }

  return $jobId;
}

sub mergeMuTect {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$rH_samplePair->{'sample'}})) {
    $jobDependency = $globalDep{$parentStep}->{$rH_samplePair->{'sample'}};
  }

  my $nbJobs = LoadConfig::getParam( $rH_cfg, 'mutect', 'approxNbJobs' );
  my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);

  my $sampleName = $rH_samplePair->{'sample'};
  my $vcfDir = 'pairedVariants/' . $sampleName."/rawMuTect/";
  my $outputDir = 'pairedVariants/' . $sampleName.'/';
  my $outputVCF = $outputDir.$sampleName.'.mutect.vcf';

  my @vcfs;
  for my $region (@{$rA_regions}) {
    push(@vcfs, $outputDir.'/rawMuTect/'.$sampleName.'.'.$region.'.mutect.vcf');
  }

  my $rO_job = VCFtools::mergeVCF($rH_cfg, \@vcfs, $outputVCF);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeMuTect", undef, 'MERGEMUTECT', $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

1;

