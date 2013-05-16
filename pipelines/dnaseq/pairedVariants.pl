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
use Breakdancer;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SnpEff;
use SubmitToCluster;
use SVtools;
use ToolShed;
use VCFtools;
use Pindel;
use Cfreec;

#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'snpAndIndelBCF'});
push(@steps, {'name' => 'mergeFilterBCF'});
push(@steps, {'name' => 'filterNStretches'});
push(@steps, {'name' => 'flagMappability'});
push(@steps, {'name' => 'snpIDAnnotation'});
push(@steps, {'name' => 'snpEffect'});
push(@steps, {'name' => 'dbNSFPAnnotation'});
push(@steps, {'name' => 'indexVCF'});
push(@steps, {'name' => 'DNAC'});
push(@steps, {'name' => 'Breakdancer'});
push(@steps, {'name' => 'Pindel'});
push(@steps, {'name' => 'ControlFreec'});

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

sub filterNStretches {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MERGEBCF_JOB_ID}';
  }

  my $sampleName = $rH_samplePair->{'sample'};
  # Use mergeFilterBCF to make sure we have the right path
  my $vcf = LoadConfig::getParam($rH_cfg, "mergeFilterBCF", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.vcf';
  my $vcfOutput = LoadConfig::getParam($rH_cfg, "filterNStretches", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.vcf';

  my $command = ToolShed::filterNStretches($rH_cfg, $sampleName, $vcf, $vcfOutput);
  my $filterNJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "filterNStretches", undef, 'FILTERN', $jobDependency, $sampleName, $command);
  return $filterNJobId;
}

sub flagMappability {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${FILTERN_JOB_ID}';
  }

  my $sampleName = $rH_samplePair->{'sample'};
  # Use mergeFilterBCF to make sure we have the right path
  my $vcf = LoadConfig::getParam($rH_cfg, "filterNStretches", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.vcf';
  my $vcfOutput = LoadConfig::getParam($rH_cfg, "flagMappability", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.vcf';

  my $command = VCFtools::annotateMappability($rH_cfg, $sampleName, $vcf, $vcfOutput);
  my $milJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "flagMappability", undef, 'MAPPABILITY', $jobDependency, $sampleName, $command);
}

sub snpIDAnnotation {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MAPPABILITY_JOB_ID}';
  }

  my $sampleName = $rH_samplePair->{'sample'};
  # Use mergeFilterBCF to make sure we have the right path
  my $vcf = LoadConfig::getParam($rH_cfg, "flagMappability", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.vcf';
  my $vcfOutput = LoadConfig::getParam($rH_cfg, "snpIDAnnotation", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.snpid.vcf';

  my $command = SnpEff::annotateDbSnp($rH_cfg, $sampleName, $vcf, $vcfOutput);
  my $snpEffJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "snpIDAnnotation", undef, 'SNPID', $jobDependency, $sampleName, $command);
}

sub snpEffect {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${SNPID_JOB_ID}';
  }

  my $sampleName = $rH_samplePair->{'sample'};
  # Use mergeFilterBCF to make sure we have the right path
  my $vcf = LoadConfig::getParam($rH_cfg, "snpIDAnnotation", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.snpid.vcf';
  my $vcfOutput = LoadConfig::getParam($rH_cfg, "snpEffect", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.snpid.snpeff.vcf';

  my $command = SnpEff::computeEffects($rH_cfg, $sampleName, $vcf, $vcfOutput);
  my $snpEffJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "snpEffect", undef, 'SNPEFF', $jobDependency, $sampleName, $command);
}

sub dbNSFPAnnotation {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${SNPEFF_JOB_ID}';
  }

  my $sampleName = $rH_samplePair->{'sample'};
  # Use mergeFilterBCF to make sure we have the right path
  my $vcf = LoadConfig::getParam($rH_cfg, "snpEffect", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.snpid.snpeff.vcf';
  my $vcfOutput = LoadConfig::getParam($rH_cfg, "dbNSFPAnnotation", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.snpid.snpeff.dbnsfp.vcf';

  my $command = SnpEff::annotateDbNSFP($rH_cfg, $sampleName, $vcf, $vcfOutput);
  my $snpEffJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "dbNSFPAnnotation", undef, 'DBNSFP', $jobDependency, $sampleName, $command);
}

sub indexVCF {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${DBNSFP_JOB_ID}';
  }

  my $sampleName = $rH_samplePair->{'sample'};
  # Use mergeFilterBCF to make sure we have the right path
  my $vcf = LoadConfig::getParam($rH_cfg, "dbNSFPAnnotation", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.snpid.snpeff.dbnsfp.vcf';
  my $vcfOutput = LoadConfig::getParam($rH_cfg, "indexVCF", 'sampleOutputRoot') . $sampleName.'/'.$sampleName.'.merged.flt.Nfilter.mil.snpid.snpeff.dbnsfp.vcf.gz';

  my $command = VCFtools::indexVCF($rH_cfg, $sampleName, $vcf, $vcfOutput);
  my $milJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "indexVCF", undef, 'INDEXVCF', $jobDependency, $sampleName, $command);
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
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.bam';
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.bam';
  my $outputDir = LoadConfig::getParam($rH_cfg, "Breakdancer", 'sampleOutputRoot') . $sampleName.'/breakdancer/';

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
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '$FILTBRD_JOB_IDS';
  }

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.bam';
  my $normalMetrics = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.sorted.dup.all.metrics.insert_size_metrics';
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.bam';
  my $tumorMetrics = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.sorted.dup.all.metrics.insert_size_metrics';
  my $outputDir = LoadConfig::getParam($rH_cfg, "Pindel", 'sampleOutputRoot') . $sampleName.'/pindel/';

  print 'mkdir -p '.$outputDir."\n";
  print "PI_JOB_IDS=\"\"\n";

  my $sampleCFGOutput = $outputDir.$sampleName.'.Pindel.conf';
  my $command = Pindel::pairedConfigFile($rH_cfg, $tumorMetrics, $normalMetrics, $tumorBam, $normalBam, $sampleCFGOutput);
  my $piCFGJobId; 
  if(defined($command) && length($command) > 0) {
    $piCFGJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "PindelCfg", undef, 'PI_CFG', $jobDependency, $sampleName, $command, LoadConfig::getParam($rH_cfg, "Pindel", 'sampleOutputRoot') . $sampleName);
    $piCFGJobId = '$'.$piCFGJobId;
  }
  my $pi2filterDep;
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $chrFile = LoadConfig::getParam($rH_cfg, "Pindel", 'referenceGenomeByChromosome') .'/chr' .$seqName.'.fa';
    my $Brdresfile = LoadConfig::getParam($rH_cfg, "Breakdancer", 'sampleOutputRoot') . $sampleName.'/breakdancer/'.$sampleName.'.brd.'.$seqName.'ctx';
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
      $piJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "pindel", $seqName, 'PI', $piCFGJobId, $sampleName, $command, LoadConfig::getParam($rH_cfg, "Pindel", 'sampleOutputRoot') . $sampleName);
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
  my $piMergeFilterJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "filterSV", $sampleName .'PI', 'FILTPI', $pi2filterDep, $sampleName, $command, LoadConfig::getParam($rH_cfg, "Pindel", 'sampleOutputRoot') . $sampleName);
  print 'FILTPI_JOB_IDS=$'.$piMergeFilterJobId."\n";
  return '$FILTPI_JOB_IDS';
}

sub ControlFreec {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rH_samplePair = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
#   if($depends > 0) {
#     $jobDependency = '$MPILEUP_JOB_IDS';
#   }

  my $sampleName = $rH_samplePair->{'sample'};
  my $normalBam = $rH_samplePair->{'normal'}.'/'.$rH_samplePair->{'normal'}.'.'.LoadConfig::getParam($rH_cfg, "ControlFreec", 'inputExtension');
  my $tumorBam = $rH_samplePair->{'tumor'}.'/'.$rH_samplePair->{'tumor'}.'.'.LoadConfig::getParam($rH_cfg, "ControlFreec", 'inputExtension');
  my $outputDir = LoadConfig::getParam($rH_cfg, "ControlFreec", 'sampleOutputRoot') . $sampleName.'/controlFREEC/';
  my $sampleConfigFile = $outputDir.'/'.$sampleName.'.freec.cfg';
  

  print 'mkdir -p '.$outputDir."\n";
  print "CONTROL_FREEC_JOB_IDS=\"\"\n";

  my $command = Cfreec::pairedFreec($rH_cfg, $tumorBam, $normalBam, $sampleConfigFile, $outputDir);
  my $cFreecJobId; 
  if(defined($command) && length($command) > 0) {
    $cFreecJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "ControlFreec", undef, 'CONTROL_FREEC', $jobDependency, $sampleName, $command, LoadConfig::getParam($rH_cfg, "ControlFreec", 'sampleOutputRoot') . $sampleName);
    $cFreecJobId = '$'.$cFreecJobId;
  }
  return $cFreecJobId;

}

1;

