#!/usr/bin/env perl

=head1 NAME

I<dnaSeq>

=head1 SYNOPSIS

  perl dnaSeq.pl -c dnaSeq.guillimin.ini -s 1 -e 21 -n project.nanuq.csv > toRun.sh

  will generate a bash script for steps 1 to 21. This script can then be executed:

  sh toRun.sh

  Options

  -c (dnaSeq.guillimin.ini) the standard configuration file for the pipeline. Templates for some cluster systems like Abacus or Guillimin may already be available at pipelines/dnaseq
  -s The start step
  -e The end step
  -n (project.nanuq.csv) the standard NANUQ read set sheet.

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

B<BVATools> Bam and Variant Analysis Tools

B<BWA> Burrows-Wheeler Aligner

B<GATK> Reads , alignment and metrics tools

B<IGVTools> igv utilities 

B<LoadConfig> Parse configuration (ini) file

B<Picard> Sort, merge tools for bam files

B<SampleSheet> Parse Nanuq sample sheet

B<SAMtools> alignment file tools

B<SequenceDictionaryParser> Parse sequence dictionnary

B<SnpEff> multiple variant annotation tools

B<SubmitToCluster> Create the submit command, control for dependencies

B<Trimmomatic> Trim and filter raw read files

B<Tools> function related to the mugqic tool shed

B<Version> Tracks app version

B<VCFtools> Manage variants vcf files

B<Metrics> Read, alignment and multiple metrics library

B<GqSeqUtils> is a library to access/launch functions from the gqSeqUtils R package


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
use Cwd 'abs_path';
use File::Basename;
use File::Path;
use Parse::Range qw(parse_range);
use POSIX;

use BVATools;
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
use Tools;
use Version;
use VCFtools;
use Metrics;
use GqSeqUtils;
#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'samToFastq', 'stepLoop' => 'sample', 'parentStep' => undef});
push(@steps, {'name' => 'trimAndAlign', 'stepLoop' => 'sample', 'parentStep' => 'samToFastq'});
push(@steps, {'name' => 'laneMetrics', 'stepLoop' => 'sample', 'parentStep' => 'trimAndAlign'});
push(@steps, {'name' => 'mergeTrimStats', 'stepLoop' => 'experiment', 'parentStep' => 'trimAndAlign'});
push(@steps, {'name' => 'symlinkRawAlignBam', 'stepLoop' => 'sample', 'parentStep' => undef});
push(@steps, {'name' => 'mergeLanes', 'stepLoop' => 'sample', 'parentStep' => 'trimAndAlign'});
push(@steps, {'name' => 'indelRealigner', 'stepLoop' => 'sample', 'parentStep' => 'mergeLanes'});
push(@steps, {'name' => 'mergeRealigned', 'stepLoop' => 'sample', 'parentStep' => 'indelRealigner'});
push(@steps, {'name' => 'fixmate', 'stepLoop' => 'sample', 'parentStep' => 'mergeRealigned'});
push(@steps, {'name' => 'markDup', 'stepLoop' => 'sample', 'parentStep' => 'fixmate'});
push(@steps, {'name' => 'recalibration', 'stepLoop' => 'sample', 'parentStep' => 'markDup'});
push(@steps, {'name' => 'metrics', 'stepLoop' => 'sample', 'parentStep' => 'recalibration'});
push(@steps, {'name' => 'callableBases', 'stepLoop' => 'sample', 'parentStep' => 'recalibration'});
push(@steps, {'name' => 'metricsLibrarySample', 'stepLoop' => 'experiment', 'parentStep' => 'metrics'});
push(@steps, {'name' => 'snpAndIndelBCF', 'stepLoop' => 'experiment', 'parentStep' => 'recalibration'});
push(@steps, {'name' => 'mergeFilterBCF', 'stepLoop' => 'experiment', 'parentStep' => 'snpAndIndelBCF'});
push(@steps, {'name' => 'filterNStretches', 'stepLoop' => 'experiment', 'parentStep' => 'mergeFilterBCF'});
push(@steps, {'name' => 'flagMappability', 'stepLoop' => 'experiment', 'parentStep' => 'filterNStretches'});
push(@steps, {'name' => 'snpIDAnnotation', 'stepLoop' => 'experiment', 'parentStep' => 'flagMappability'});
push(@steps, {'name' => 'snpEffect', 'stepLoop' => 'experiment', 'parentStep' => 'snpIDAnnotation'});
push(@steps, {'name' => 'dbNSFPAnnotation', 'stepLoop' => 'experiment', 'parentStep' => 'snpEffect'});
push(@steps, {'name' => 'metricsSNV', 'stepLoop' => 'experiment', 'parentStep' => 'snpEffect'});
push(@steps, {'name' => 'deliverable' , 'stepLoop' => 'experiment' , 'parentStep' => ['mergeTrimStats','metricsLibrarySample','metricsSNV']});
push(@steps, {'name' => 'fullPileup', 'stepLoop' => 'sample', 'parentStep' => 'recalibration'});

#--------------------
# PODS
#--------------------
## Here starts the pipeline steps documentation, please change it accordingly any time you add/remove/modify a step

=head1 RNASEQ PIPELINE STEPS

The standard variant discovery pipeline performs the following steps:

B<trimAndAlign> :  Raw reads quality trimming and removing of Illumina adapters is performed using trimmomatic. The filtered reads are aligned to a reference genome. The alignment is done per lane of sequencing. The alignment software used is bwa, two different algorithms are available: bwa aln (backtrack) and bwa mem.

B<laneMetrics> :  Aligned reads are duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score) will be marked as a duplicate in the .bam file. Marking duplicates is done with picard software. Metrics per lane are produced at this step.

B<mergeTrimStats> : The trim statistics per lane/sample are merged at this step.

B<mergeLanes> : Bam files per sample are merged in one file. Merge is done using the Picard software.

B<indelRealigner> : Insertion and deletion realignment is performed on regions where multiple base mismatches are preferred over indels by the aligner since it can appear to be less costly by the algorithm. Such regions will introduce false positive variant calls which may be filtered out by realigning those regions properly. Realignment is done using the GATK software. The reference genome is divided by a number regions given by the nbRealignJobs parameter.

B<mergeRealigned> : Merging the regions of realigned reads per sample is done using Picard.

B<fixmate> : Fixing the read mates. Once local regions are realigned, the read mate coordinates of the aligned reads need to be recalculated since the reads are realigned at positions that differ from their original alignment. Fixing the read mate positions is done using the Picard software.

B<markDup> : Marking duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score) will be marked as a duplicate in the .bam file. Marking duplicates is done using the Picard software.

B<recalibration> : Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration, the quality scores in the QUAL field in each read in the output BAM are more accurate in that the reported quality score is closer to its actual probability of mismatching the reference genome. Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle and sequence context, and by doing so provides not only more accurate quality scores but also more widely dispersed ones.

B<metrics> : Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads, Median, mean and standard deviation of insert sizes of reads after alignment, mean coverage over exons (mean number of reads per base position), percentage of bases covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads) whole genome
percentage of bases covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads) for specific targets (CCDS regions are used in human samples). A TDF (.tdf) coverage track is also generated at this step for easy visualization of coverage in the IGV browser.

B<metricsLibrarySample> : Merge metrics. Read metrics per sample are merged at this step.

B<snpAndIndelBCF> : Mpileup and Variant calling. Variants (SNPs and INDELs) are called using samtools mpileup. bcftools view is used to produce binary bcf files.  

B<mergeFilterBCF> : bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step.  The output of bcftools is fed to varfilter, which does an additional filtering of the variants and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls for all samples in the experiment. 

B<filterNStretches> : The final .vcf files are filtered for long 'N' INDELs which are sometimes introduced and cause excessive memory usage by downstream tools.

B<flagMappability> : Mappability annotation. An in-house database identifies regions in which reads are confidently mapped to the reference genome.

B<snpIDAnnotation>:  dbSNP annotation. The .vcf files are annotated for dbSNP using the software SnpSift (from the SNPEff suite).

B<snpEffect>: Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software. SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

B<dbNSFPAnnotation> Additional SVN annotations. Provides extra information about SVN by using numerous published databases. Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information) and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy) and other function annotations) .

B<metricsSNV> : Metrics SNV. Multiple metrics associated to annotations and effect prediction are generated at this step: change rate by chromosome, changes by type, effects by impact, effects by functional class, counts by effect, counts by genomic region, SNV quality, coverage, InDel lengths, base changes,  transition-transversion rates, summary of allele frequencies, codon changes, amino acid changes, changes per chromosome, change rates.

B<deliverable> : Generating the standard report. A summary html report is automatically generated by the pipeline. This report contains description of the sequencing experiment as well as a detailed presentation of the pipeline steps and results. Various Quality Control (QC) summary statistics are included in the report and additional QC analysis is accessible for download directly through the report. The report includes also the main references of the software and methods used during the analysis, together with the full list of parameters passed to the pipeline main script.

B<fullPileup>: Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format. One packaged mpileup file is created per sample/chromosome.

=cut end of documentation


my %globalDep;
for my $stepName (@steps) {
  $globalDep{$stepName -> {'name'} } ={};
}


# Global scope variables
my $configFile;
my $workDirectory = getcwd();


&main();

sub printUsage {
  print "Version: ".$Version::version."\n";
  print "\nUsage: perl ".$0." -c config.ini -s start -e end -n SampleSheet.csv\n";
  print "\t-c  config file\n";
  print "\t-s  step range e.g. '1,3', '2-5', '1,4-7,10'\n";
  print "\t-n  nanuq sample sheet\n";
  #print "\t-d  First dependency (optional)\n";
  print "\n";
  print "Steps:\n";
  for(my $idx=0; $idx < @steps; $idx++) {
    print "".($idx+1).'- '.$steps[$idx]->{'name'}."\n";
  }
  print "\n";
}

sub main {
  my %opts;
  getopts('c:s:n:d:', \%opts);
  
  if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'n'})) {
    printUsage();
    exit(1);
  }

  my $firstDependency = $opts{'d'};
  my %cfg = LoadConfig->readConfigFile($opts{'c'});
  my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($opts{'n'});
  my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);
  $configFile =  abs_path($opts{'c'});

  my @sampleNames = keys %{$rHoAoH_sampleInfo};

  print STDERR "Samples: ".scalar(@sampleNames)."\n";
  if(defined($firstDependency) && length($firstDependency) > 0) {
    for my $sampleName (@sampleNames) {
      $globalDep{"default"}->{$sampleName} = $firstDependency;
    }
  }

  # List user-defined step index range.
  # Shift 1st position to 0 instead of 1
  my @stepRange = map($_ - 1, parse_range($opts{'s'}));

  SubmitToCluster::initPipeline;

  my $currentStep;
  my $lastStepId = $stepRange[$#stepRange];
  for(my $idx=0; $idx < @sampleNames; $idx++){
    my $sampleName = $sampleNames[$idx];
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};

    for $currentStep (@stepRange) {
      my $fname = $steps[$currentStep]->{'name'};
      my $subref = \&$fname;

      if ($steps[$currentStep]->{'stepLoop'} eq 'sample') {
        # Tests for the first step in the list. Used for dependencies.
        my $jobIdVar = &$subref($currentStep, \%cfg, $sampleName, $rAoH_sampleLanes, $rAoH_seqDictionary); 
        $globalDep{$fname}->{$sampleName} = $jobIdVar;

        if(defined($jobIdVar) && $currentStep == $lastStepId) {
          print "FINAL_STEP_".$idx.'_JOB_IDS='.$jobIdVar."\n";
        }
      }
    }
  }

  for $currentStep (@stepRange) {
    if($steps[$currentStep]->{'stepLoop'} eq 'experiment') {
      my $fname = $steps[$currentStep]->{'name'};
      my $subref = \&$fname;

      my $jobIdVar = &$subref($currentStep, \%cfg, $rHoAoH_sampleInfo, $rAoH_seqDictionary);
      $globalDep{$fname}->{'experiment'} = $jobIdVar;
    }
  }
  
  
  my $jobId = "";
  if($steps[$lastStepId]->{'stepLoop'} eq 'experiment') {
    if(defined($globalDep{$steps[$lastStepId]->{'name'}}->{'experiment'})) {
      $jobId = $globalDep{$steps[$lastStepId]->{'name'}}->{'experiment'};
    }
  }
  else {
    my @finalIds;
    for(my $idx=0; $idx < @sampleNames; $idx++){
      push(@finalIds, '${FINAL_STEP_'.$idx.'_JOB_IDS}');
    }
    $jobId = join(':', @finalIds);
  }
  print 'export FINAL_STEP_JOB_IDS='.$jobId."\n";
  
}

sub samToFastq {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = 'default';
  if(defined($parentStep) && defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  print "SAMTOFASTQ_JOB_IDS=\"\"\n";
  my $setJobId = 0;
  my $first = 1;

  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    # Check if bam is defined
    $rH_laneInfo->{'bam'} or die "Error in dnaSeq::samToFastq: BAM file is not defined for sample run lane $sampleName " . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "!";

    my $rO_job;
    my $baseDirectory = "$sampleName/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rawDirectory = LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir', 1, 'dirpath') . "/" . $baseDirectory;
    my $input1 = $rawDirectory . "/" . $rH_laneInfo->{'bam'};

    if ($input1) {
      if ($rH_laneInfo->{'runType'} eq "SINGLE_END") {
        my $outputFastq1 = $input1;
        $outputFastq1 =~ s/\.bam$/.single.fastq.gz/;
        $rO_job = Picard::samToFastq($rH_cfg, $input1, $outputFastq1);
        $rH_laneInfo->{'read1File'} = basename($outputFastq1);
      } elsif ($rH_laneInfo->{'runType'} eq "PAIRED_END") {
        my $outputFastq1 = $input1;
        my $outputFastq2 = $input1;
        $outputFastq1 =~ s/\.bam$/.pair1.fastq.gz/;
        $outputFastq2 =~ s/\.bam$/.pair2.fastq.gz/;
        $rO_job = Picard::samToFastq($rH_cfg, $input1, $outputFastq1, $outputFastq2);
        $rH_laneInfo->{'read1File'} = basename($outputFastq1);
        $rH_laneInfo->{'read2File'} = basename($outputFastq2);
      } else {
        die "Error in dnaSeq::samToFastq: unknown run type (can be 'SINGLE_END' or 'PAIRED_END' only): " . $rH_laneInfo->{'runType'};
      }
    }
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "samToFastq", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, "SAMTOFASTQ", $jobDependency, $sampleName, $rO_job);
      if ($first == 1) {
        print "SAMTOFASTQ_JOB_IDS=" . $rO_job->getCommandJobId(0) . "\n";
        $first = 0;
      } else {
        print "SAMTOFASTQ_JOB_IDS=\${SAMTOFASTQ_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0) . "\n";
      }
      $setJobId = 1;
    }
  }

  if ($setJobId) {
    return "\${SAMTOFASTQ_JOB_IDS}";
  } else {
    return undef;
  }
}

sub trimAndAlign {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($parentStep) && defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  print "BWA_JOB_IDS=\"\"\n";
  my $setJobId = 0;
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $rgId = $rH_laneInfo->{'libraryBarcode'} . "_" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rgSampleName = $rH_laneInfo->{'name'};
    my $rgLibrary = $rH_laneInfo->{'libraryBarcode'};
    my $rgPlatformUnit = $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rgCenter = undef;

    my $outputDir = 'reads/'.$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    print 'mkdir -p '.$outputDir."\n";
    my $ro_trimJob = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo, $outputDir);
    SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM', $jobDependency, $sampleName, $ro_trimJob);
    
    my $trimDependency = $ro_trimJob->getCommandJobId(0);
    if(defined($jobDependency) && (!defined($trimDependency) || length ($trimDependency) == 0)) {
      $trimDependency = $jobDependency;
    }

    my $outputAlnDir = 'alignment/'.$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    print 'mkdir -p '.$outputAlnDir."\n";
    my $outputAlnPrefix = $outputAlnDir.'/'.$sampleName.'.'.$rH_laneInfo->{'libraryBarcode'};

    my $useMem = LoadConfig::getParam($rH_cfg, 'aln', 'aligner') eq 'mem';
    if(!$useMem) {
      $rgCenter = LoadConfig::getParam( $rH_cfg, 'aln', 'bwaInstitution' );
      my $ro_bwaJob = BWA::aln($rH_cfg, $sampleName, $ro_trimJob->getOutputFileHash()->{PAIR1_OUTPUT}, $ro_trimJob->getOutputFileHash()->{PAIR2_OUTPUT},$ro_trimJob->getOutputFileHash()->{SINGLE1_OUTPUT}, $outputAlnPrefix, $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter);
      if(!$ro_bwaJob->isUp2Date()) {
        if($ro_bwaJob->getNbCommands() == 3) {
            SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read1.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ1ALN', $trimDependency, $sampleName, $ro_bwaJob, 0 );
            SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read2.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ2ALN', $trimDependency, $sampleName, $ro_bwaJob, 1 );
            SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'sampe.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA', $ro_bwaJob->getCommandJobId(0).LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$ro_bwaJob->getCommandJobId(1), $sampleName, $ro_bwaJob, 2 );
            print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$ro_bwaJob->getCommandJobId(2)."\n";
            $setJobId = 1;
        }
        else {
          SubmitToCluster::printSubmitCmd($rH_cfg, "aln", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READALN', $trimDependency, $sampleName, $ro_bwaJob, 0 );
          SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'samse.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA',  $ro_bwaJob->getCommandJobId(1), $sampleName, $ro_bwaJob, 1 );
          print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$ro_bwaJob->getCommandJobId(1)."\n";
          $setJobId = 1;
        } 
      }
    }
    else {
      $rgCenter = LoadConfig::getParam( $rH_cfg, 'mem', 'bwaInstitution' );
      my $ro_bwaJob = BWA::mem($rH_cfg, $sampleName, $ro_trimJob->getOutputFileHash()->{PAIR1_OUTPUT}, $ro_trimJob->getOutputFileHash()->{PAIR2_OUTPUT},$ro_trimJob->getOutputFileHash()->{SINGLE1_OUTPUT}, $outputAlnPrefix, $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter);
      if(!$ro_bwaJob->isUp2Date()) {
        SubmitToCluster::printSubmitCmd($rH_cfg, "mem", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA_MEM', $trimDependency, $sampleName, $ro_bwaJob);
        print 'BWA_JOB_IDS='.$ro_bwaJob->getCommandJobId(0)."\n";
        $setJobId = 1;
      }
    }
  }

  if($setJobId ==0){
    return undef;
  }

  return '$BWA_JOB_IDS';
}

sub laneMetrics {
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

  print "LANE_METRICS_JOB_IDS=\"\"\n";
  my $first=1;
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $directory = 'alignment/'.$sampleName."/run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'}."/";
    my $sortedLaneBamFile = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.bam';
    my $sortedLaneDupBamFile = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.dup.bam';
    my $outputMetrics = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.dup.metrics';

    my $rO_job = Picard::markDup($rH_cfg, $sampleName, $sortedLaneBamFile, $sortedLaneDupBamFile, $outputMetrics);
    if(!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'LANEMARKDUP', $jobDependency, $sampleName, $rO_job);
      if($first == 1) {
        print 'LANE_METRICS_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
        $first = 0;
      }
      else {
        print 'LANE_METRICS_JOB_IDS=${LANE_METRICS_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'laneMarkDup', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
      }
    }

    $outputMetrics = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.dup.metrics';
    my $rO_collectMetricsJob = Picard::collectMetrics($rH_cfg, $sortedLaneBamFile, $outputMetrics);
    if(!$rO_collectMetricsJob->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "collectMetrics", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'COLLECTMETRICS', $jobDependency, $sampleName, $rO_collectMetricsJob);
      if($first == 1) {
        print 'LANE_METRICS_JOB_IDS='.$rO_collectMetricsJob->getCommandJobId(0)."\n";
      }
      else {
        print 'LANE_METRICS_JOB_IDS=${LANE_METRICS_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'collectMetrics', 'clusterDependencySep').$rO_collectMetricsJob->getCommandJobId(0)."\n";
      }
    }

    $outputMetrics = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.dup.metrics.nodup.targetCoverage.txt';
    my $coverageBED = BVATools::resolveSampleBED($rH_cfg, $rH_laneInfo);
    my $rO_coverageJob = BVATools::depthOfCoverage($rH_cfg, $sortedLaneBamFile, $outputMetrics, $coverageBED);
    if(!$rO_coverageJob->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "depthOfCoverage", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'LANEDEPTHOFCOVERAGE', $jobDependency, $sampleName, $rO_coverageJob);
      if($first == 1) {
        print 'LANE_METRICS_JOB_IDS='.$rO_coverageJob->getCommandJobId(0)."\n";
      }
      else {
        print 'LANE_METRICS_JOB_IDS=${LANE_METRICS_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'clusterDependencySep').$rO_coverageJob->getCommandJobId(0)."\n";
      }
    }
  }

  return '${LANE_METRICS_JOB_IDS}';
}

sub mergeTrimStats {
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
  my $rO_job = Metrics::mergeTrimmomaticStats($rH_cfg,  $libraryType, $pattern, $folder, $ouputFile);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "trimMetrics", undef, 'TRIMMETRICS', $trimmingDependency, undef, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}


sub symlinkRawAlignBam {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my @inputBams;

  

  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $baseDirectory = "$sampleName/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rawBam = abs_path(LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir', 1, 'dirpath')) . "/" . $baseDirectory . "/" . $rH_laneInfo->{'bam'};

    my $alignmentDirectory = 'alignment/' . $baseDirectory . "/";
    my $alignedBam = $alignmentDirectory . "/" . $rH_laneInfo->{'name'} . '.' . $rH_laneInfo->{'libraryBarcode'} . '.sorted.bam';

    mkpath $alignmentDirectory;

    if (-l $alignedBam) {
      warn "[Warning] Symbolic link $alignedBam already exists! Skipping.\n";
    } elsif (-e $rawBam and symlink($rawBam, $alignedBam)) {
      warn "Created symbolic link $alignedBam successfully.\n";
    } else {
      die "[Error] Can't create symbolic link $alignedBam to target $rawBam!\n";
    }
  }

  return undef;
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

  my @inputBams;
  my $outputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam';
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $directory = 'alignment/'.$sampleName."/run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'}."/";
    my $sortedLaneBamFile = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.bam';

    push(@inputBams, $sortedLaneBamFile);
  }

  my $rO_job = Picard::mergeFiles($rH_cfg, $sampleName, \@inputBams, $outputBAM);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeLanes", undef, 'MERGELANES', $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
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

  my $nbRealignJobs = LoadConfig::getParam($rH_cfg, 'indelRealigner', 'nbRealignJobs', 1, 'int');
  if($nbRealignJobs > 50) {
    warn "Number of realign jobs is >50. This is usually much. Anything beyond 20 can be problematic.\n";
  }


  my $jobId;
  if($nbRealignJobs <= 1) {
    my $rO_job = GATK::realign($rH_cfg, $sampleName, 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam', undef, 'alignment/'.$sampleName.'/realign/all');
    if(!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", undef, 'REALIGN', $jobDependency, $sampleName, $rO_job);
      print 'REALIGN_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
      $jobId = '${REALIGN_JOB_IDS}';
    }
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
    my $firstJob = 1;
    for my $seqName (@chrToProcess) {
      push(@excludeList, $seqName);
      my $rO_job = GATK::realign($rH_cfg, $sampleName, 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam', $seqName, 'alignment/'.$sampleName.'/realign/'.$seqName, $processUnmapped);
      if($processUnmapped == 1) {
        $processUnmapped = 0;
      }

      if(!$rO_job->isUp2Date()) {
        SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", $seqName, 'REALIGN', $jobDependency, $sampleName, $rO_job);
        if($firstJob) {
          $jobId = $rO_job->getCommandJobId(0);
          print 'REALIGN_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
          $firstJob = 0;
        }
        else {
          print 'REALIGN_JOB_IDS=${REALIGN_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
        }
        $jobId = '${REALIGN_JOB_IDS}';
      }
    }

    my $rO_job = GATK::realign($rH_cfg, $sampleName, 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.bam', undef, 'alignment/'.$sampleName.'/realign/others', 1, \@excludeList);
    if(!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", 'others', 'REALIGN', $jobDependency, $sampleName, $rO_job);
      print 'REALIGN_JOB_IDS=${REALIGN_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
      $jobId = '${REALIGN_JOB_IDS}';
    }
  }
  return $jobId;
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

  my $jobId;
  my @inputBams;
  my $outputBAM = 'alignment/'.$sampleName.'/'.$sampleName.'.realigned.qsorted.bam';

  my $nbRealignJobs = LoadConfig::getParam($rH_cfg, 'indelRealigner', 'nbRealignJobs', 1, 'int');
  my $rO_job;
  if($nbRealignJobs <= 1) {
    my $command = 'if [ ! -e '.$outputBAM.' ]; then ln -s alignment/'.$sampleName.'/realign/all.bam '.$outputBAM.'; fi;';
    $rO_job = new Job(0);
    $rO_job->addCommand($command);
  }
  else {
    #Keep space for the exclude realignment at the end.
    $nbRealignJobs--;
    for (my $idx=0; $idx < $nbRealignJobs; $idx++) {
      my $input = 'alignment/'.$sampleName.'/realign/'.$rAoH_seqDictionary->[$idx]->{'name'}.'.bam';
      push(@inputBams, $input);
    }
    push(@inputBams, 'alignment/'.$sampleName.'/realign/others.bam');

    $rO_job = Picard::mergeFiles($rH_cfg, $sampleName, \@inputBams, $outputBAM);
  }

  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeRealign", undef, 'MERGEREALIGN', $jobDependency, $sampleName, $rO_job);
    $jobId = $rO_job->getCommandJobId(0);
  }

  return $jobId;
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

  my $rO_job = Picard::fixmate($rH_cfg, $inputBAM, $outputBAM);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "fixmate", undef, 'FIXMATE', $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
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

  my $rO_job = Picard::markDup($rH_cfg, $sampleName, $inputBAM, $outputBAM, $outputMetrics);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'MARKDUP', $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
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

  my $rO_job = GATK::recalibration($rH_cfg, $sampleName, $inputBAM, $outputBAM);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "recalibration", undef, 'RECAL', $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
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

  my $bamFile = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.bam';
  my $jobId;

  # Collect metrics
  my $outputMetrics = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.all.metrics';
  my $rO_collectMetricsJob = Picard::collectMetrics($rH_cfg, $bamFile, $outputMetrics);
  if(!$rO_collectMetricsJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "collectMetrics", undef, 'COLLECTMETRICS', $jobDependency, $sampleName, $rO_collectMetricsJob);
    if(!defined($jobId)) {
      $jobId='$METRICS_JOBS';
      print 'METRICS_JOBS='.$rO_collectMetricsJob->getCommandJobId(0)."\n";
    }
  }
  
  # Compute genome coverage
  my $outputPrefix = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.all.coverage';
  my $rO_genomeCoverageJob = GATK::genomeCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  if(!$rO_genomeCoverageJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "genomeCoverage", undef, 'GENOMECOVERAGE', $jobDependency, $sampleName, $rO_genomeCoverageJob);
    if(!defined($jobId)) {
      $jobId='$METRICS_JOBS';
      print 'METRICS_JOBS='.$rO_genomeCoverageJob->getCommandJobId(0)."\n";;
    }
    else {
      print 'METRICS_JOBS=${METRICS_JOBS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_genomeCoverageJob->getCommandJobId(0)."\n";
    }
  }

  my $output = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.coverage.tsv';
  my $coverageBED = BVATools::resolveSampleBED($rH_cfg, $rAoH_sampleLanes->[0]);
  my $rO_coverageJob = BVATools::depthOfCoverage($rH_cfg, $bamFile, $output, $coverageBED);
  if(!$rO_coverageJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "depthOfCoverage", undef, 'SAMPLE_COVERAGE', $jobDependency, $sampleName, $rO_coverageJob);
    if(!defined($jobId)) {
      $jobId='$METRICS_JOBS';
      print 'METRICS_JOBS='.$rO_coverageJob->getCommandJobId(0)."\n";;
    }
    else {
      print 'METRICS_JOBS=${METRICS_JOBS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_coverageJob->getCommandJobId(0)."\n";
    }
  }
  # Compute CCDS coverage
  $outputPrefix = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.CCDS.coverage';
  my $rO_targetCoverageJob = GATK::targetCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  if(!$rO_targetCoverageJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "targetCoverage", undef, 'TARGETCOVERAGE', $jobDependency, $sampleName, $rO_targetCoverageJob);
    if(!defined($jobId)) {
      $jobId='$METRICS_JOBS';
      print 'METRICS_JOBS='.$rO_targetCoverageJob->getCommandJobId(0)."\n";;
    }
    else {
      print 'METRICS_JOBS=${METRICS_JOBS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_targetCoverageJob->getCommandJobId(0)."\n";
    }
  }

  # Generate IGV track
  my $rO_igvtoolsTDFJob = IGVTools::computeTDF($rH_cfg, $bamFile);
  if(!$rO_igvtoolsTDFJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "computeTDF", undef, 'IGVTOOLS', $jobDependency, $sampleName, $rO_igvtoolsTDFJob);
    if(!defined($jobId)) {
      $jobId='$METRICS_JOBS';
      print 'METRICS_JOBS='.$rO_igvtoolsTDFJob->getCommandJobId(0)."\n";
    }
    else {
      print 'METRICS_JOBS=${METRICS_JOBS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_igvtoolsTDFJob->getCommandJobId(0)."\n";
    }
  }

  # Compute flags
  $output = 'alignment/'.$sampleName.'/'.$sampleName.'.sorted.dup.recal.bam.flagstat';
  my $rO_flagstatJob = SAMtools::flagstat($rH_cfg, $bamFile, $output);
  if(!$rO_flagstatJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "flagstat", undef, 'FLAGSTAT', $jobDependency, $sampleName, $rO_flagstatJob);
    if(!defined($jobId)) {
      $jobId='$METRICS_JOBS';
      print 'METRICS_JOBS='.$rO_flagstatJob->getCommandJobId(0)."\n";
    }
    else {
      print 'METRICS_JOBS=${METRICS_JOBS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_flagstatJob->getCommandJobId(0)."\n";
    }
  }

  return $jobId;
}

sub callableBases {
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
  my $jobId=undef;

  my $outputPrefix = 'alignment/'.$sampleName.'/'.$sampleName;
  my $rO_job = GATK::callableBases($rH_cfg, $bamFile, $outputPrefix);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "callableBases", undef, 'CALLABLE_BASES', $jobDependency, $sampleName, $rO_job);
    $jobId='$CALLABLE_BASES_JOB='.$rO_job->getCommandJobId(0)."\n";
  }

  return $jobId;
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

  my $rO_job = Metrics::mergeSampleDnaStats($rH_cfg,  $experimentType, $folder, $ouputFile);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "sampleMetrics", undef, 'SAMPLEMETRICS', $metricsDependency, undef, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
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
  my $jobId;
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $outputPerSeq = $outputDir.$sampleName.'.'.$seqName.'.mpileup.gz';
    my $rO_job = SAMtools::rawmpileup($rH_cfg, $sampleName, $bamFile, $seqName, $outputPerSeq);
    if(!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "rawmpileup", $seqName, 'RAW_MPILEUP', $jobDependency, $sampleName, $rO_job);
      if(!defined($jobId)) {
        $jobId = '${RAW_MPILEUP_JOB_IDS}';
        print 'RAW_MPILEUP_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
      }
      else {
        print 'RAW_MPILEUP_JOB_IDS=${RAW_MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
      }
    }

    $catCommand .= $outputPerSeq .' ';
  }

  if(defined($jobId)) {
    my $output = $outputDir.$sampleName.'.mpileup.gz';
    $catCommand .= '| gzip -c --best > '.$output;

    my $rO_job = new Job(0);
    $rO_job->addCommand($catCommand);

    SubmitToCluster::printSubmitCmd($rH_cfg, "rawmpileup_cat", undef, 'RAW_MPILEUP_CAT', '${RAW_MPILEUP_JOB_IDS}', $sampleName, $rO_job);
    return $rO_job->getCommandJobId(0);
  }
  return undef;
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

  my $jobId;
  my $nbJobs = LoadConfig::getParam($rH_cfg, 'mpileup', 'approxNbJobs', 0, 'int');
  if (defined($nbJobs) && $nbJobs > 1) {
    my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);
    for my $region (@{$rA_regions}) {
      my $rO_job = SAMtools::mpileup($rH_cfg, 'allSamples', \@inputFiles, $region, $outputDir);
      if(!$rO_job->isUp2Date()) {
        $region =~ s/:/_/g;
        SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $region, 'MPILEUP', $jobDependencies, 'allSamples', $rO_job);
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
    my $region = undef;
    my $rO_job = SAMtools::mpileup($rH_cfg, 'allSamples', \@inputFiles, $region, $outputDir);
    if(!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $region, 'MPILEUP', $jobDependencies, 'allSamples', $rO_job);
      if(!defined($jobId)) {
        $jobId = '${MPILEUP_JOB_IDS}';
        print 'MPILEUP_JOB_IDS='.$rO_job->getCommandJobId(0)."\n";
      }
      else {
        print 'MPILEUP_JOB_IDS=${MPILEUP_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$rO_job->getCommandJobId(0)."\n";
      }
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

  my $nbJobs = LoadConfig::getParam($rH_cfg, 'mpileup', 'approxNbJobs');
  my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);

  my $bcfDir = 'variants/rawBCF/';
  my $outputDir = 'variants/';

  my $rO_job = SAMtools::mergeFilterBCF($rH_cfg, 'allSamples', $bcfDir, $outputDir, $rA_regions);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeFilterBCF", undef, 'MERGEBCF', $jobDependency, 'allSamples', $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub filterNStretches {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $inputVCF = 'variants/allSamples.merged.flt.vcf';
  my $outputVCF = 'variants/allSamples.merged.flt.NFiltered.vcf';

  my $rO_job = Tools::filterNStretches($rH_cfg, 'allSamples', $inputVCF, $outputVCF);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "filterNStretches", undef, 'FILTERN', $jobDependency, 'allSamples', $rO_job);
  }
  return $rO_job->getCommandJobId(0);
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

  my $inputVCF = 'variants/allSamples.merged.flt.NFiltered.vcf';
  my $outputVCF = 'variants/allSamples.merged.flt.mil.vcf';
  my $rO_job = VCFtools::annotateMappability($rH_cfg, $inputVCF, $outputVCF);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "flagMappability", undef, 'MAPPABILITY', $jobDependency, 'allSamples', $rO_job);
  }
  return $rO_job->getCommandJobId(0);
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

  my $rO_job = SnpEff::annotateDbSnp($rH_cfg, $inputVCF, $vcfOutput);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "snpIDAnnotation", undef, 'SNPID', $jobDependency, 'allSamples', $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

sub snpEffect {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $inputVCF = 'variants/allSamples.merged.flt.mil.snpId.vcf';
  my $vcfOutput = 'variants/allSamples.merged.flt.mil.snpId.snpeff.vcf';

  my $rO_job = SnpEff::computeEffects($rH_cfg, $inputVCF, $vcfOutput, 1);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "snpEffect", undef, 'SNPEFF', $jobDependency, 'allSamples', $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

sub dbNSFPAnnotation {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }
  

  my $inputVCF = 'variants/allSamples.merged.flt.mil.snpId.snpeff.vcf';
  my $vcfOutput = 'variants/allSamples.merged.flt.mil.snpId.snpeff.dbnsfp.vcf';
  
  my $rO_job = SnpEff::annotateDbNSFP($rH_cfg, $inputVCF, $vcfOutput);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "dbNSFPAnnotation", undef, 'DBNSFP', $jobDependency, 'allSamples', $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

sub metricsSNV {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $vcfDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};

  if(defined($globalDep{$parentStep}->{'experiment'})){
      $vcfDependency .= $globalDep{$parentStep}->{'experiment'};
  }
  
  
  my $inputVCF = 'variants/allSamples.merged.flt.mil.snpId.vcf';
  my $outputFile = 'metrics/allSamples.merged.flt.mil.snpId.snpeff.vcf.part.changeRate.tsv';
  my $listFiles='variants/allSamples.merged.flt.mil.snpId.snpeff.vcf.statsFile.txt';

  my $rO_job_changeRate = Metrics::svnStatsChangeRate($rH_cfg, $inputVCF, $outputFile, $listFiles);
  if(!$rO_job_changeRate->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "metricsSNV", undef, 'CHANGERATE', $vcfDependency , 'allSamples', $rO_job_changeRate);
  }
  
  my $outputBaseName='metrics/allSamples.SNV';
  my $rO_job_Graph = Metrics::svnStatsGetGraph($rH_cfg, $listFiles,$outputBaseName);
  if(!$rO_job_Graph->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "metricsSNV", "2", 'SNV_GRAPH', $rO_job_changeRate->getCommandJobId(0) , 'allSamples', $rO_job_Graph);
  }
  return $rO_job_Graph->getCommandJobId(0);
}


sub deliverable {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;


  my $reportDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};

  my $jobDependencies = "";
  foreach my $stepName (@{$parentStep}) {
    if(defined($globalDep{$stepName}->{'experiment'})){
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$globalDep{$stepName}->{'experiment'};
    }
  }
  if(length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }
  $reportDependency = $jobDependencies;


  my $rO_job = GqSeqUtils::clientReport($rH_cfg,  $configFile, $workDirectory, 'DNAseq') ;

  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "deliverable", 'REPORT', 'DNAREPORT', $reportDependency , 'allSamples', $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}


1;
