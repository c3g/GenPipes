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

B<BVATools> BAM and Variant Analysis Tools

B<BWA> Burrows-Wheeler Aligner

B<GATK> Reads, alignment and metrics tools

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
use strict qw(vars subs);
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
use File::Path 'make_path';
use POSIX;

use BVATools;
use BWA;
use GATK;
use GqSeqUtils;
use IGVTools;
use LoadConfig;
use Metrics;
use Picard;
use Pipeline;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SnpEff;
use SubmitToCluster;
use Trimmomatic;
use Tools;
use Version;
use VCFtools;
#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'bamToFastq', 'loop' => 'readSet', 'parentSteps' => []});
push(@steps, {'name' => 'trim', 'loop' => 'readSet', 'parentSteps' => ['bamToFastq']});
push(@steps, {'name' => 'bwaAln1', 'loop' => 'readSet', 'parentSteps' => ['trim']});
push(@steps, {'name' => 'bwaAln2', 'loop' => 'readSet', 'parentSteps' => ['trim']});
push(@steps, {'name' => 'bwaSortSam', 'loop' => 'readSet', 'parentSteps' => ['bwaAln1', 'bwaAln2']});
push(@steps, {'name' => 'bwaMem', 'loop' => 'readSet', 'parentSteps' => ['trim']});
push(@steps, {'name' => 'laneMetrics', 'loop' => 'sample', 'parentSteps' => ['bwaSortSam', 'bwaMem']});
push(@steps, {'name' => 'mergeTrimStats', 'loop' => 'global', 'parentSteps' => ['bwaSortSam', 'bwaMem']});
push(@steps, {'name' => 'mergeLanes', 'loop' => 'sample', 'parentSteps' => ['bwaSortSam', 'bwaMem']});
push(@steps, {'name' => 'indelRealigner', 'loop' => 'sample', 'parentSteps' => ['mergeLanes']});
push(@steps, {'name' => 'mergeRealigned', 'loop' => 'sample', 'parentSteps' => ['indelRealigner']});
push(@steps, {'name' => 'fixmate', 'loop' => 'sample', 'parentSteps' => ['mergeRealigned']});
push(@steps, {'name' => 'markDup', 'loop' => 'sample', 'parentSteps' => ['fixmate']});
push(@steps, {'name' => 'recalibration', 'loop' => 'sample', 'parentSteps' => ['markDup']});
push(@steps, {'name' => 'metrics', 'loop' => 'sample', 'parentSteps' => ['recalibration']});
push(@steps, {'name' => 'metricsLibrarySample', 'loop' => 'global', 'parentSteps' => ['metrics']});
push(@steps, {'name' => 'snpAndIndelBCF', 'loop' => 'global', 'parentSteps' => ['recalibration']});
push(@steps, {'name' => 'mergeFilterBCF', 'loop' => 'global', 'parentSteps' => ['snpAndIndelBCF']});
push(@steps, {'name' => 'filterNStretches', 'loop' => 'global', 'parentSteps' => ['mergeFilterBCF']});
push(@steps, {'name' => 'flagMappability', 'loop' => 'global', 'parentSteps' => ['filterNStretches']});
push(@steps, {'name' => 'snpIDAnnotation', 'loop' => 'global', 'parentSteps' => ['flagMappability']});
push(@steps, {'name' => 'snpEffect', 'loop' => 'global', 'parentSteps' => ['snpIDAnnotation']});
push(@steps, {'name' => 'dbNSFPAnnotation', 'loop' => 'global', 'parentSteps' => ['snpEffect']});
push(@steps, {'name' => 'metricsSNV', 'loop' => 'global', 'parentSteps' => ['snpEffect']});
push(@steps, {'name' => 'deliverable', 'loop' => 'global', 'parentSteps' => ['mergeTrimStats', 'metricsLibrarySample', 'metricsSNV']});
push(@steps, {'name' => 'fullPileup', 'loop' => 'sample', 'parentSteps' => ['recalibration']});

#--------------------
# PODS
#--------------------
## Here starts the pipeline steps documentation, please change it accordingly any time you add/remove/modify a step

=head1 DNASEQ PIPELINE STEPS

The standard variant discovery pipeline performs the following steps:

B<trimAndAlign> : Raw reads quality trimming and removing of Illumina adapters is performed using trimmomatic. The filtered reads are aligned to a reference genome. The alignment is done per lane of sequencing. The alignment software used is bwa, two different algorithms are available: bwa aln (backtrack) and bwa mem.

B<laneMetrics> : Aligned reads are duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score) will be marked as a duplicate in the .bam file. Marking duplicates is done with picard software. Metrics per lane are produced at this step.

B<mergeTrimStats> : The trim statistics per lane/sample are merged at this step.

B<mergeLanes> : BAM files per sample are merged in one file. Merge is done using the Picard software.

B<indelRealigner> : Insertion and deletion realignment is performed on regions where multiple base mismatches are preferred over indels by the aligner since it can appear to be less costly by the algorithm. Such regions will introduce false positive variant calls which may be filtered out by realigning those regions properly. Realignment is done using the GATK software. The reference genome is divided by a number regions given by the nbRealignJobs parameter.

B<mergeRealigned> : Merging the regions of realigned reads per sample is done using Picard.

B<fixmate> : Fixing the read mates. Once local regions are realigned, the read mate coordinates of the aligned reads need to be recalculated since the reads are realigned at positions that differ from their original alignment. Fixing the read mate positions is done using the Picard software.

B<markDup> : Marking duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score) will be marked as a duplicate in the .bam file. Marking duplicates is done using the Picard software.

B<recalibration> : Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration, the quality scores in the QUAL field in each read in the output BAM are more accurate in that the reported quality score is closer to its actual probability of mismatching the reference genome. Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle and sequence context, and by doing so provides not only more accurate quality scores but also more widely dispersed ones.

B<metrics> : Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads, Median, mean and standard deviation of insert sizes of reads after alignment, mean coverage over exons (mean number of reads per base position), percentage of bases covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads) whole genome
percentage of bases covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads) for specific targets (CCDS regions are used in human samples). A TDF (.tdf) coverage track is also generated at this step for easy visualization of coverage in the IGV browser.

B<metricsLibrarySample> : Merge metrics. Read metrics per sample are merged at this step.

B<snpAndIndelBCF> : Mpileup and Variant calling. Variants (SNPs and INDELs) are called using samtools mpileup. bcftools view is used to produce binary bcf files.

B<mergeFilterBCF> : bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step. The output of bcftools is fed to varfilter, which does an additional filtering of the variants and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls for all samples in the experiment.

B<filterNStretches> : The final .vcf files are filtered for long 'N' INDELs which are sometimes introduced and cause excessive memory usage by downstream tools.

B<flagMappability> : Mappability annotation. An in-house database identifies regions in which reads are confidently mapped to the reference genome.

B<snpIDAnnotation>: dbSNP annotation. The .vcf files are annotated for dbSNP using the software SnpSift (from the SNPEff suite).

B<snpEffect>: Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software. SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

B<dbNSFPAnnotation> Additional SVN annotations. Provides extra information about SVN by using numerous published databases. Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information) and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy) and other function annotations) .

B<metricsSNV> : Metrics SNV. Multiple metrics associated to annotations and effect prediction are generated at this step: change rate by chromosome, changes by type, effects by impact, effects by functional class, counts by effect, counts by genomic region, SNV quality, coverage, InDel lengths, base changes, transition-transversion rates, summary of allele frequencies, codon changes, amino acid changes, changes per chromosome, change rates.

B<deliverable> : Generating the standard report. A summary html report is automatically generated by the pipeline. This report contains description of the sequencing experiment as well as a detailed presentation of the pipeline steps and results. Various Quality Control (QC) summary statistics are included in the report and additional QC analysis is accessible for download directly through the report. The report includes also the main references of the software and methods used during the analysis, together with the full list of parameters passed to the pipeline main script.

B<fullPileup>: Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format. One packaged mpileup file is created per sample/chromosome.

=cut end of documentation


my %globalDep;
for my $stepName (@steps) {
  $globalDep{$stepName -> {'name'} } = {};
}


# Global scope variables
my $configFile;
my $workDirectory = getcwd();


&main();

sub debug {
  my $message = shift;

  my $debug = 1;    # Set to 1 to display debug messages, 0 to keep output silent

  $debug and print STDERR "[DEBUG] $message\n";
}

sub getUsage {
  my $usage = <<END;
MUGQIC Pipeline DNA-Seq Version: $Version::version

Usage: perl $0 -c CONFIG_FILE -r STEP_RANGE -s SAMPLE_FILE [-o OUTPUT_DIR]
  -c  config file
  -r  step range e.g. "1-5", "3,6,7", "2,4-8"
  -s  sample file
  -o  output directory (default: current)

Steps:
END

  # List and number step names
  for (my $i = 1; $i <= @steps; $i++) {
    $usage .= $i . "- " . $steps[$i - 1]->{'name'} . "\n";
  }

  return $usage;
}

sub main {
  my %opts;
  getopts('c:r:s:o:', \%opts);

  if (!defined($opts{'c'}) || !defined($opts{'r'}) || !defined($opts{'s'})) {
    die getUsage();
  }

  my $stepRange = $opts{'r'};
  my $outputDirectory = $opts{'o'};
  my %cfg = LoadConfig->readConfigFile($opts{'c'});
  my $sampleFile = $opts{'s'};
  my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);

  my $pipeline = Pipeline->new(\@steps, $sampleFile, $outputDirectory);

  # Go through steps and create global or readSet jobs accordingly
  foreach my $step ($pipeline->getStepsByRange($stepRange)) {
    my $stepName = $step->getName();
    debug "main: processing step $stepName";

    # ReadSet step creates 1 job per readSet per sample
    if ($step->getLoop() eq 'readSet') {
      foreach my $sample (@{$pipeline->getSamples()}) {
        foreach my $readSet (@{$sample->getReadSets()}) {
          debug "main: processing read set " . $readSet->getName();

          my $rO_job = &$stepName(\%cfg, $readSet, $rAoH_seqDictionary);
          if ($rO_job) {
            $rO_job->setLoopTags([$sample->getName(), $readSet->getName()]);
            debug "main: readSet job " . join(".", @{$rO_job->getLoopTags()});
            #debug "main: readSet job input files " . join(" ", @{$rO_job->getInputFiles()});
            $step->submitJob(\%cfg, $rO_job);
          }
        }
      }
    # Sample step creates 1 job per sample
    } elsif ($step->getLoop() eq 'sample') {
      foreach my $sample (@{$pipeline->getSamples()}) {
        debug "main: processing sample " . $sample->getName();

        my $rO_job = &$stepName(\%cfg, $sample, $rAoH_seqDictionary);
        if ($rO_job) {
          $rO_job->setLoopTags([$sample->getName()]);
          debug "main: sample job " . join(".", @{$rO_job->getLoopTags()});
          $step->submitJob(\%cfg, $rO_job);
        }
      }
    # Global step creates 1 job only
    } else {
      my $rO_job = &$stepName(\%cfg, $rAoH_seqDictionary);
      if ($rO_job) {
        $rO_job->setLoopTags([]);
        $step->submitJob(\%cfg, $rO_job);
      }
    }
  }
}

sub bamToFastq {
  my $rH_cfg = shift;
  my $readSet = shift;
  my $rAoH_seqDictionary = shift;

  my $rO_job;

  if ($readSet->getBAM() and not($readSet->getFASTQ1())) {
    if ($readSet->getRunType() eq "PAIRED_END") {
      my $FASTQ1 = $readSet->getBAM();
      $FASTQ1 =~ s/\.bam$/.pair1.fastq.gz/;
      $readSet->setFASTQ1($FASTQ1);

      my $FASTQ2 = $readSet->getBAM();
      $FASTQ2 =~ s/\.bam$/.pair2.fastq.gz/;
      $readSet->setFASTQ2($FASTQ2);

      $rO_job = Picard::samToFastq($rH_cfg, $readSet->getBAM(), $readSet->getFASTQ1(), $readSet->getFASTQ2());
    } elsif ($readSet->getRunType() eq "SINGLE_END") {
      my $FASTQ1 = $readSet->getBAM();
      $FASTQ1 =~ s/\.bam$/.single.fastq.gz/;
      $readSet->setFASTQ1($FASTQ1);

      $rO_job = Picard::samToFastq($rH_cfg, $readSet->getBAM(), $readSet->getFASTQ1());
    }
  }
  return $rO_job;
}

sub trim {
  my $rH_cfg = shift;
  my $readSet = shift;

  my $trimFilePrefix = "\$WORK_DIR/trim/" . $readSet->getSample()->getName() . "/" . $readSet->getName() . ".trim.";

  my $pairedOutput1;
  my $unpairedOutput1;
  my $pairedOutput2;
  my $unpairedOutput2;
  my $singleOutput;

  if ($readSet->getRunType() eq "PAIRED_END") {
    $pairedOutput1 = $trimFilePrefix . "pair1.fastq.gz";
    $unpairedOutput1 = $trimFilePrefix . "single1.fastq.gz";
    $pairedOutput2 = $trimFilePrefix . "pair2.fastq.gz";
    $unpairedOutput2 = $trimFilePrefix . "single2.fastq.gz";
  } elsif ($readSet->getRunType() eq "SINGLE_END") {
    $singleOutput = $trimFilePrefix . "single.fastq.gz";
  } else {
    die "[Error] In trim, unknown runType: \"" . $readSet->getRunType() . "\"!";
  }

  return Trimmomatic::trim(
    $rH_cfg,
    $readSet->getFASTQ1(),
    $readSet->getFASTQ2(),
    $pairedOutput1,
    $unpairedOutput1,
    $pairedOutput2,
    $unpairedOutput2,
    $singleOutput,
    $readSet->getQualityOffset(),
    $trimFilePrefix . "out",
    $trimFilePrefix . "stats.csv"
  );
}

sub bwaAln1 {
  my $rH_cfg = shift;
  my $readSet = shift;

  my $inDbFasta = LoadConfig::getParam($rH_cfg, 'aln', 'bwaRefIndex', 1, 'filepath');
  my $inQueryFastq;
  my $trimFilePrefix = "\$WORK_DIR/trim/" . $readSet->getSample()->getName() . "/" . $readSet->getName() . ".trim.";
  my $outSai;
  my $saiFilePrefix = "\$WORK_DIR/alignment/" . $readSet->getSample()->getName() . "/" . $readSet->getName();

  if ($readSet->getRunType() eq "PAIRED_END") {
    $inQueryFastq = $trimFilePrefix . "pair1.fastq.gz";
    $outSai = $saiFilePrefix . ".pair1.sai";
  } elsif ($readSet->getRunType() eq "SINGLE_END") {
    $inQueryFastq = $trimFilePrefix . "single.fastq.gz";
    $outSai = $saiFilePrefix . ".single.sai";
  } else {
    die "[Error] In bwaAln1, unknown runType: \"" . $readSet->getRunType() . "\"!";
  }

  return BWA::aln(
    $rH_cfg,
    $inDbFasta,
    $inQueryFastq,
    $outSai
  );
}

sub bwaAln2 {
  my $rH_cfg = shift;
  my $readSet = shift;

  my $inDbFasta = LoadConfig::getParam($rH_cfg, 'aln', 'bwaRefIndex', 1, 'filepath');
  my $inQueryFastq;
  my $trimFilePrefix = "\$WORK_DIR/trim/" . $readSet->getSample()->getName() . "/" . $readSet->getName() . ".trim.";
  my $outSai;
  my $saiFilePrefix = "\$WORK_DIR/alignment/" . $readSet->getSample()->getName() . "/" . $readSet->getName();

  if ($readSet->getRunType() eq "PAIRED_END") {
    $inQueryFastq = $trimFilePrefix . "pair2.fastq.gz";
    $outSai = $saiFilePrefix . ".pair2.sai";

    return BWA::aln(
      $rH_cfg,
      $inDbFasta,
      $inQueryFastq,
      $outSai
    );
  } elsif ($readSet->getRunType() eq "SINGLE_END") {
    # No second alignment for SINGLE_END reads
    return undef;
  } else {
    die "[Error] In bwaAln2, unknown runType: \"" . $readSet->getRunType() . "\"!";
  }
}


sub bwaSortSam {
  my $rH_cfg = shift;
  my $readSet = shift;

  my $inDbFasta = LoadConfig::getParam($rH_cfg, 'aln', 'bwaRefIndex', 1, 'filepath');

  my $trimFilePrefix = "\$WORK_DIR/trim/" . $readSet->getSample()->getName() . "/" . $readSet->getName() . ".trim.";
  my $alignmentFilePrefix = "\$WORK_DIR/alignment/" . $readSet->getSample()->getName() . "/" . $readSet->getName();

  my $rgLibrary = $readSet->getLibrary();
  my $rgPlatformUnit = $readSet->getRun() . "_" . $readSet->getLane();
  my $rgId = $rgLibrary . "_" . $rgPlatformUnit;
  my $rgSampleName = $readSet->getSample()->getName();
  my $rgCenter = LoadConfig::getParam($rH_cfg, 'aln', 'bwaInstitution');
  my $readGroup = "'" . '@RG\tID:' . $rgId . '\tSM:' . $rgSampleName . '\tLB:' . $rgLibrary . '\tPU:run' . $rgPlatformUnit . '\tCN:' . $rgCenter . '\tPL:Illumina' . "'";


  my $rO_bwaSamJob;
  my $rO_sortSamJob;

  if ($readSet->getRunType() eq "PAIRED_END") {
    $rO_bwaSamJob = BWA::sampe(
      $rH_cfg,
      $inDbFasta,
      $alignmentFilePrefix . ".pair1.sai",
      $alignmentFilePrefix . ".pair2.sai",
      $trimFilePrefix . "pair1.fastq.gz",
      $trimFilePrefix . "pair2.fastq.gz",
      undef,
      $readGroup
    );
  } elsif ($readSet->getRunType() eq "SINGLE_END") {
    $rO_bwaSamJob = BWA::samse(
      $rH_cfg,
      $inDbFasta,
      $alignmentFilePrefix . ".single.sai",
      $trimFilePrefix . "single.fastq.gz",
      undef,
      $readGroup
    );
  } else {
    die "[Error] In bwaSortSam, unknown runType: \"" . $readSet->getRunType() . "\"!";
  }

  $rO_sortSamJob = Picard::sortSam(
    $rH_cfg,
    "/dev/stdin",
    $alignmentFilePrefix . ".sorted.bam",
    "coordinate"
  );

  return Job::pipe($rO_bwaSamJob, $rO_sortSamJob);
}

sub bwaTmp {
  my $rH_cfg = shift;
  my $readSet = shift;


  my $rgLibrary = $readSet->getLibrary();
  my $rgPlatformUnit = $readSet->getRun() . "_" . $readSet->getLane();
  my $rgId = $rgLibrary . "_" . $rgPlatformUnit;
  my $rgSampleName = $readSet->getSample()->getName();
  my $rgCenter = undef;


  my $outputAlnDir = "alignment/" . $readSet->getSample()->getName() . "/" . $readSet->getName();
  print "mkdir -p " . $outputAlnDir . "\n";
  my $outputAlnPrefix = $outputAlnDir . "/" . $readSet->getName();



  my $rH_laneInfo = shift;
  my $trimDir;
  my $trimDependency;
  my $sampleName;
  my $setJobId;
  my $bwaPair1Input = undef;
  my $bwaPair2Input = undef;
  my $bwaSingleInput = undef;

  if ($rH_laneInfo->{'read2File'}) {
    $bwaPair1Input = $trimDir . "/" . basename($rH_laneInfo->{'read1File'});
    $bwaPair1Input =~ s/\.pair1\.fastq\.gz$/.trim.pair1.fastq.gz/;
    $bwaPair2Input = $trimDir . "/" . basename($rH_laneInfo->{'read2File'});
    $bwaPair2Input =~ s/\.pair2\.fastq\.gz$/.trim.pair2.fastq.gz/;
  } else {
    $bwaSingleInput = $trimDir . "/" . basename($rH_laneInfo->{'read1File'});
    $bwaSingleInput =~ s/\.single\.fastq\.gz$/.trim.single.fastq.gz/;
  }

  my $aligner = LoadConfig::getParam($rH_cfg, 'aln', 'aligner');
  if ($aligner eq "aln") {
    $rgCenter = LoadConfig::getParam($rH_cfg, 'aln', 'bwaInstitution');
    my $rO_bwaJob = BWA::aln(
      $rH_cfg,
      $sampleName,
      $bwaPair1Input,
      $bwaPair2Input,
      $bwaSingleInput,
      $outputAlnPrefix,
      $rgId,
      $rgSampleName,
      $rgLibrary,
      $rgPlatformUnit,
      $rgCenter
    );
    if (!$rO_bwaJob->isUp2Date()) {
      if ($rO_bwaJob->getNbCommands() == 3) {
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", "read1." . $rgPlatformUnit, "READ1ALN", $trimDependency, $sampleName, $rO_bwaJob, 0);
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", "read2." . $rgPlatformUnit, "READ2ALN", $trimDependency, $sampleName, $rO_bwaJob, 1);
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", "sampe." . $rgPlatformUnit, "BWA", $rO_bwaJob->getCommandJobId(0) . LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep') . $rO_bwaJob->getCommandJobId(1), $sampleName, $rO_bwaJob, 2);
        print "BWA_JOB_IDS=\${BWA_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep') . $rO_bwaJob->getCommandJobId(2) . "\n";
          $setJobId = 1;
      } else {
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", $rgPlatformUnit, "READALN", $trimDependency, $sampleName, $rO_bwaJob, 0);
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", "samse." . $rgPlatformUnit, "BWA", $rO_bwaJob->getCommandJobId(1), $sampleName, $rO_bwaJob, 1);
        print "BWA_JOB_IDS=\${BWA_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_bwaJob->getCommandJobId(1) . "\n";
      }
    }
  } else {
    $rgCenter = LoadConfig::getParam($rH_cfg, 'mem', 'bwaInstitution');
    my $rO_bwaJob = BWA::mem(
      $rH_cfg,
      $sampleName,
      $bwaPair1Input,
      $bwaPair2Input,
      $bwaSingleInput,
      $outputAlnPrefix,
      $rgId,
      $rgSampleName,
      $rgLibrary,
      $rgPlatformUnit,
      $rgCenter);
    if (!$rO_bwaJob->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "mem", $rgPlatformUnit, "BWA_MEM", $trimDependency, $sampleName, $rO_bwaJob);
      print "BWA_JOB_IDS=" . $rO_bwaJob->getCommandJobId(0) . "\n";
      $setJobId = 1;
    }
  }
}

sub laneMetrics {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  print "LANE_METRICS_JOB_IDS=\"\"\n";
  my $first = 1;
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $directory = "alignment/" . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    my $sortedLaneBAMFile = $directory . $rH_laneInfo->{'name'} . "." . $rH_laneInfo->{'libraryBarcode'} . '.sorted.bam';
    my $sortedLaneDupBAMFile = $directory . $rH_laneInfo->{'name'} . "." . $rH_laneInfo->{'libraryBarcode'} . '.sorted.dup.bam';
    my $outputMetrics = $directory . $rH_laneInfo->{'name'} . "." . $rH_laneInfo->{'libraryBarcode'} . '.sorted.dup.metrics';

    my $rO_job = Picard::markDup($rH_cfg, $sampleName, $sortedLaneBAMFile, $sortedLaneDupBAMFile, $outputMetrics);
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, "LANEMARKDUP", $jobDependency, $sampleName, $rO_job);
      if ($first == 1) {
        print "LANE_METRICS_JOB_IDS=" . $rO_job->getCommandJobId(0) . "\n";
        $first = 0;
      } else {
        print "LANE_METRICS_JOB_IDS=\${LANE_METRICS_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'laneMarkDup', 'clusterDependencySep') . $rO_job->getCommandJobId(0) . "\n";
      }
    }

    $outputMetrics = $directory . $rH_laneInfo->{'name'} . "." . $rH_laneInfo->{'libraryBarcode'} . ".sorted.dup.metrics";
    my $rO_collectMetricsJob = Picard::collectMetrics($rH_cfg, $sortedLaneBAMFile, $outputMetrics);
    if (!$rO_collectMetricsJob->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "collectMetrics", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, "COLLECTMETRICS", $jobDependency, $sampleName, $rO_collectMetricsJob);
      if ($first == 1) {
        print "LANE_METRICS_JOB_IDS=" . $rO_collectMetricsJob->getCommandJobId(0) . "\n";
      } else {
        print "LANE_METRICS_JOB_IDS=\${LANE_METRICS_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'collectMetrics', 'clusterDependencySep') . $rO_collectMetricsJob->getCommandJobId(0) . "\n";
      }
    }

    $outputMetrics = $directory . $rH_laneInfo->{'name'} . "." . $rH_laneInfo->{'libraryBarcode'} . ".sorted.dup.metrics.nodup.targetCoverage.txt";
    my $coverageBED = BVATools::resolveSampleBED($rH_cfg, $rH_laneInfo);
    my $rO_coverageJob = BVATools::depthOfCoverage($rH_cfg, $sortedLaneBAMFile, $outputMetrics, $coverageBED);
    if (!$rO_coverageJob->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "depthOfCoverage", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, "LANEDEPTHOFCOVERAGE", $jobDependency, $sampleName, $rO_coverageJob);
      if ($first == 1) {
        print "LANE_METRICS_JOB_IDS=" . $rO_coverageJob->getCommandJobId(0) . "\n";
      } else {
        print "LANE_METRICS_JOB_IDS=\${LANE_METRICS_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'clusterDependencySep') . $rO_coverageJob->getCommandJobId(0) . "\n";
      }
    }
  }

  return "\${LANE_METRICS_JOB_IDS}";
}

sub mergeTrimStats {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rAoH_seqDictionary = shift;

  my $libraryType = undef;
  my $fkey = (keys %{$rHoAoH_sampleInfo})[0];
  my @fvals = @{$rHoAoH_sampleInfo->{$fkey}};
  my $finfo = $fvals[0];
  if ($finfo->{'runType'} eq "SINGLE_END") {
    $libraryType = 'single';
  } elsif ($finfo->{'runType'} eq "PAIRED_END") {
    $libraryType = 'paired';
  }
  my $trimmingDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};

  my @sampleNames = keys %{$rHoAoH_sampleInfo};
  my $jobDependencies = "";
  for (my $idx = 0; $idx < @sampleNames; $idx++) {
    my $sampleName = $sampleNames[$idx];
    if (defined($globalDep{$parentStep}->{$sampleName})) {
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $globalDep{$parentStep}->{$sampleName};
    }
  }
  if (length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }
  $trimmingDependency = $jobDependencies;

  my $folder = "reads";
  my $pattern = "trim.stats.csv";
  my $ouputFile = "metrics/trimming.stats";
  print "mkdir -p metrics\n";
  my $rO_job = Metrics::mergeTrimmomaticStats($rH_cfg, $libraryType, $pattern, $folder, $ouputFile);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "trimMetrics", undef, "TRIMMETRICS", $trimmingDependency, undef, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}


sub mergeLanes {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my @inputBAMs;
  my $outputBAM = "alignment/$sampleName/$sampleName.sorted.bam";
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $baseDirectory = "$sampleName/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rawDirectory = LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir', 1, 'dirpath') . "/" . $baseDirectory;
    my $inputBAM = abs_path($rawDirectory . "/" . $rH_laneInfo->{'read1File'});
    my $alignDirectory = "alignment/$baseDirectory";
    my $sortedLaneBAMFile = $alignDirectory . "/" . $rH_laneInfo->{'name'} . "." . $rH_laneInfo->{'libraryBarcode'} . ".sorted.bam";

    if ($inputBAM =~ /\.bam$/) {
      make_path $alignDirectory;

      if (-l $sortedLaneBAMFile) {
        warn "[Warning] Symbolic link $sortedLaneBAMFile already exists! Skipping.\n";
      } elsif (-f $inputBAM and symlink($inputBAM, $sortedLaneBAMFile)) {
        warn "Created symbolic link $sortedLaneBAMFile successfully.\n";
      } else {
        die "[Error] Can't create symbolic link $sortedLaneBAMFile to target $inputBAM!\n";
      }
    }

    push(@inputBAMs, $sortedLaneBAMFile);
  }

  my $rO_job = Picard::mergeFiles($rH_cfg, $sampleName, \@inputBAMs, $outputBAM);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeLanes", undef, "MERGELANES", $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub indelRealigner {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }
  print "mkdir -p alignment/" . $sampleName . "/realign\n";

  my $nbRealignJobs = LoadConfig::getParam($rH_cfg, 'indelRealigner', 'nbRealignJobs', 1, 'int');
  if ($nbRealignJobs > 50) {
    warn "Number of realign jobs is >50. This is usually much. Anything beyond 20 can be problematic.\n";
  }


  my $jobId;
  if ($nbRealignJobs <= 1) {
    my $rO_job = GATK::realign($rH_cfg, $sampleName, "alignment/" . $sampleName . "/" . $sampleName . ".sorted.bam", undef, "alignment/" . $sampleName . "/realign/all");
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", undef, "REALIGN", $jobDependency, $sampleName, $rO_job);
      print "REALIGN_JOB_IDS=" . $rO_job->getCommandJobId(0) . "\n";
      $jobId = "\${REALIGN_JOB_IDS}";
    }
  } else {
    #Keep space for the exclude realignment at the end.
    $nbRealignJobs--;
    my @chrToProcess;
    for (my $idx = 0; $idx < $nbRealignJobs; $idx++) {
      push(@chrToProcess, $rAoH_seqDictionary->[$idx]->{'name'});
    }

    print "REALIGN_JOB_IDS=\"\"\n";
    my $processUnmapped = 1;
    my @excludeList;
    my $firstJob = 1;
    for my $seqName (@chrToProcess) {
      push(@excludeList, $seqName);
      my $rO_job = GATK::realign($rH_cfg, $sampleName, "alignment/" . $sampleName . "/" . $sampleName . ".sorted.bam", $seqName, "alignment/" . $sampleName . "/realign/" . $seqName, $processUnmapped);
      if ($processUnmapped == 1) {
        $processUnmapped = 0;
      }

      if (!$rO_job->isUp2Date()) {
        SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", $seqName, "REALIGN", $jobDependency, $sampleName, $rO_job);
        if ($firstJob) {
          $jobId = $rO_job->getCommandJobId(0);
          print "REALIGN_JOB_IDS=" . $rO_job->getCommandJobId(0) . "\n";
          $firstJob = 0;
        } else {
          print "REALIGN_JOB_IDS=\${REALIGN_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0) . "\n";
        }
        $jobId = "\${REALIGN_JOB_IDS}";
      }
    }

    my $rO_job = GATK::realign($rH_cfg, $sampleName, "alignment/" . $sampleName . "/" . $sampleName . ".sorted.bam", undef, "alignment/" . $sampleName . "/realign/others", 1, \@excludeList);
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", "others", "REALIGN", $jobDependency, $sampleName, $rO_job);
      print "REALIGN_JOB_IDS=\${REALIGN_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0) . "\n";
      $jobId = "\${REALIGN_JOB_IDS}";
    }
  }
  return $jobId;
}

sub mergeRealigned {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $jobId;
  my @inputBAMs;
  my $outputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".realigned.qsorted.bam";

  my $nbRealignJobs = LoadConfig::getParam($rH_cfg, 'indelRealigner', 'nbRealignJobs', 1, 'int');
  my $rO_job;
  if ($nbRealignJobs <= 1) {
    my $command = 'if [ ! -e ' . $outputBAM . ' ]; then ln -s alignment/' . $sampleName . '/realign/all.bam ' . $outputBAM . '; fi;';
    $rO_job = new Job(0);
    $rO_job->addCommand($command);
  } else {
    #Keep space for the exclude realignment at the end.
    $nbRealignJobs--;
    for (my $idx = 0; $idx < $nbRealignJobs; $idx++) {
      my $input = "alignment/" . $sampleName . "/realign/" . $rAoH_seqDictionary->[$idx]->{'name'} . ".bam";
      push(@inputBAMs, $input);
    }
    push(@inputBAMs, "alignment/" . $sampleName . "/realign/others.bam");

    $rO_job = Picard::mergeFiles($rH_cfg, $sampleName, \@inputBAMs, $outputBAM);
  }

  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeRealign", undef, "MERGEREALIGN", $jobDependency, $sampleName, $rO_job);
    $jobId = $rO_job->getCommandJobId(0);
  }

  return $jobId;
}

sub fixmate {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".realigned.qsorted.bam";
  my $outputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".matefixed.sorted.bam";

  my $rO_job = Picard::fixmate($rH_cfg, $inputBAM, $outputBAM);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "fixmate", undef, "FIXMATE", $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub markDup {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".matefixed.sorted.bam";
  my $outputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.bam";
  my $outputMetrics = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.metrics";

  my $rO_job = Picard::markDup($rH_cfg, $sampleName, $inputBAM, $outputBAM, $outputMetrics);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, "MARKDUP", $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub recalibration {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.bam";
  my $outputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup";

  my $rO_job = GATK::recalibration($rH_cfg, $sampleName, $inputBAM, $outputBAM);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "recalibration", undef, "RECAL", $jobDependency, $sampleName, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub metrics {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $bamFile = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.recal.bam";
  my $jobId;

  # Collect metrics
  my $outputMetrics = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.recal.all.metrics";
  my $rO_collectMetricsJob = Picard::collectMetrics($rH_cfg, $bamFile, $outputMetrics);
  if (!$rO_collectMetricsJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "collectMetrics", undef, "COLLECTMETRICS", $jobDependency, $sampleName, $rO_collectMetricsJob);
    if (!defined($jobId)) {
      $jobId = "\$METRICS_JOBS";
      print "METRICS_JOBS=" . $rO_collectMetricsJob->getCommandJobId(0) . "\n";
    }
  }

  # Compute genome coverage
  my $outputPrefix = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.recal.all.coverage";
  my $rO_genomeCoverageJob = GATK::genomeCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  if (!$rO_genomeCoverageJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "genomeCoverage", undef, "GENOMECOVERAGE", $jobDependency, $sampleName, $rO_genomeCoverageJob);
    if (!defined($jobId)) {
      $jobId = "\$METRICS_JOBS";
      print "METRICS_JOBS=" . $rO_genomeCoverageJob->getCommandJobId(0) . "\n";;
    } else {
      print "METRICS_JOBS=\${METRICS_JOBS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_genomeCoverageJob->getCommandJobId(0) . "\n";
    }
  }

  # Compute CCDS coverage
  $outputPrefix = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.recal.CCDS.coverage";
  my $rO_targetCoverageJob = GATK::targetCoverage($rH_cfg, $sampleName, $bamFile, $outputPrefix);
  if (!$rO_targetCoverageJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "targetCoverage", undef, "TARGETCOVERAGE", $jobDependency, $sampleName, $rO_targetCoverageJob);
    if (!defined($jobId)) {
      $jobId = "\$METRICS_JOBS";
      print "METRICS_JOBS=" . $rO_targetCoverageJob->getCommandJobId(0) . "\n";;
    } else {
      print "METRICS_JOBS=\${METRICS_JOBS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_targetCoverageJob->getCommandJobId(0) . "\n";
    }
  }

  # Generate IGV track
  my $rO_igvtoolsTDFJob = IGVTools::computeTDF($rH_cfg, $bamFile);
  if (!$rO_igvtoolsTDFJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "computeTDF", undef, "IGVTOOLS", $jobDependency, $sampleName, $rO_igvtoolsTDFJob);
    if (!defined($jobId)) {
      $jobId = "\$METRICS_JOBS";
      print "METRICS_JOBS=" . $rO_igvtoolsTDFJob->getCommandJobId(0) . "\n";
    } else {
      print "METRICS_JOBS=\${METRICS_JOBS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_igvtoolsTDFJob->getCommandJobId(0) . "\n";
    }
  }

  # Compute flags
  my $output = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.recal.bam.flagstat";
  my $rO_flagstatJob = SAMtools::flagstat($rH_cfg, $bamFile, $output);
  if (!$rO_flagstatJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "flagstat", undef, "FLAGSTAT", $jobDependency, $sampleName, $rO_flagstatJob);
    if (!defined($jobId)) {
      $jobId = "\$METRICS_JOBS";
      print "METRICS_JOBS=" . $rO_flagstatJob->getCommandJobId(0) . "\n";
    } else {
      print "METRICS_JOBS=\${METRICS_JOBS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_flagstatJob->getCommandJobId(0) . "\n";
    }
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
  for (my $idx = 0; $idx < @sampleNames; $idx++) {
    my $sampleName = $sampleNames[$idx];
    if (defined($globalDep{$parentStep}->{$sampleName})) {
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $globalDep{$parentStep}->{$sampleName};
    }
  }
  if (length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }
  $metricsDependency = $jobDependencies;

  my $folder = "alignment/";
  my $ouputFile = "metrics/SampleMetrics.stats";
  my $experimentType = LoadConfig::getParam($rH_cfg, 'default', 'experimentType');
  print "mkdir -p metrics\n";

  my $rO_job = Metrics::mergeSampleDnaStats($rH_cfg, $experimentType, $folder, $ouputFile);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "sampleMetrics", undef, "SAMPLEMETRICS", $metricsDependency, undef, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}


#push(@steps, {'name' => 'countTelomere'});
sub fullPileup {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if (defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $bamFile = "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.recal.bam";
  my $outputDir = "alignment/" . $sampleName . "/mpileup/";

  print "mkdir -p " . $outputDir . "\n";
  print "RAW_MPILEUP_JOB_IDS=\"\"\n";
  my $catCommand = "zcat ";
  my $jobId;
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $outputPerSeq = $outputDir . $sampleName . "." . $seqName . ".mpileup.gz";
    my $rO_job = SAMtools::rawmpileup($rH_cfg, $sampleName, $bamFile, $seqName, $outputPerSeq);
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "rawmpileup", $seqName, "RAW_MPILEUP", $jobDependency, $sampleName, $rO_job);
      if (!defined($jobId)) {
        $jobId = "\${RAW_MPILEUP_JOB_IDS}";
        print "RAW_MPILEUP_JOB_IDS=" . $rO_job->getCommandJobId(0) . "\n";
      } else {
        print "RAW_MPILEUP_JOB_IDS=\${RAW_MPILEUP_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0) . "\n";
      }
    }

    $catCommand .= $outputPerSeq . " ";
  }

  if (defined($jobId)) {
    my $output = $outputDir . $sampleName . ".mpileup.gz";
    $catCommand .= "| gzip -c --best > " . $output;

    my $rO_job = new Job(0);
    $rO_job->addCommand($catCommand);

    SubmitToCluster::printSubmitCmd($rH_cfg, "rawmpileup_cat", undef, "RAW_MPILEUP_CAT", "\${RAW_MPILEUP_JOB_IDS}", $sampleName, $rO_job);
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
  for (my $idx = 0; $idx < @sampleNames; $idx++) {
    my $sampleName = $sampleNames[$idx];
    if (defined($globalDep{$parentStep}->{$sampleName})) {
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $globalDep{$parentStep}->{$sampleName};
    }
    push(@inputFiles, "alignment/" . $sampleName . "/" . $sampleName . ".sorted.dup.recal.bam");
  }
  if (length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }

  my $outputDir = "variants/rawBCF/";
  print "mkdir -p " . $outputDir . "\n";
  print "MPILEUP_JOB_IDS=\"\"\n";

  my $jobId;
  my $nbJobs = LoadConfig::getParam($rH_cfg, 'mpileup', 'approxNbJobs', 0, 'int');
  if (defined($nbJobs) && $nbJobs > 1) {
    my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);
    for my $region (@{$rA_regions}) {
      my $rO_job = SAMtools::mpileup($rH_cfg, "allSamples", \@inputFiles, $region, $outputDir);
      if (!$rO_job->isUp2Date()) {
        $region =~ s/:/_/g;
        SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $region, "MPILEUP", $jobDependencies, "allSamples", $rO_job);
        if (!defined($jobId)) {
          $jobId = "\${MPILEUP_JOB_IDS}";
          print "MPILEUP_JOB_IDS=" . $rO_job->getCommandJobId(0) . "\n";
        } else {
          print "MPILEUP_JOB_IDS=\${MPILEUP_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0) . "\n";
        }
      }
    }
  } else {
    my $region = undef;
    my $rO_job = SAMtools::mpileup($rH_cfg, "allSamples", \@inputFiles, $region, $outputDir);
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "mpileup", $region, "MPILEUP", $jobDependencies, "allSamples", $rO_job);
      if (!defined($jobId)) {
        $jobId = "\${MPILEUP_JOB_IDS}";
        print "MPILEUP_JOB_IDS=" . $rO_job->getCommandJobId(0) . "\n";
      } else {
        print "MPILEUP_JOB_IDS=\${MPILEUP_JOB_IDS}" . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0) . "\n";
      }
    }
  }

  return $jobId;
}

sub generateApproximateWindows {
  my $nbJobs = shift;
  my $rAoH_seqDictionary = shift;

  my @retVal;
  if ($nbJobs <= scalar(@{$rAoH_seqDictionary})) {
    for my $rH_seqInfo (@$rAoH_seqDictionary) {
      push(@retVal, $rH_seqInfo->{'name'} . ":1-" . $rH_seqInfo->{'size'});
    }
  } else {
    $nbJobs -= @$rAoH_seqDictionary;
    my $totalSize = 0;
    for my $rH_seqInfo (@$rAoH_seqDictionary) {
      $totalSize += $rH_seqInfo->{'size'};
    }
    my $approxWindow = floor($totalSize / $nbJobs);

    for my $rH_seqInfo (@$rAoH_seqDictionary) {
      for (my $idx = 1; $idx <= $rH_seqInfo->{'size'}; $idx += $approxWindow) {
        my $end = $idx+$approxWindow-1;
        if ($end > $rH_seqInfo->{'size'}) {
          $end = $rH_seqInfo->{'size'};
        }

        my $region = $rH_seqInfo->{'name'} . ":" . $idx . "-" . $end;
        push(@retVal, $region);
      }
    }
  }

#  my @retVal;
#  $totalSize = 0;
#  my $currentregion = "";
#  my $approxWindowRemaining = $approxWindow;
#  for my $rH_seqInfo (@$rAoH_seqDictionary) {
#    for (my $idx = 1; $idx <= $rH_seqInfo->{'size'}; $idx += $approxWindowRemaining) {
#      my $end = $idx+$snvWindow-1;
#      my $hitEnd = 0;
#      if ($end > $rH_seqInfo->{'size'}) {
#        $end = $rH_seqInfo->{'size'};
#        $hitEnd = 1;
#        $approxWindowRemaining -= ($end - $idx) +1;
#      }
#
#      my $region = $rH_seqInfo->{'name'} . ":" . $idx . "-" . $end;
#      if (length($currentregion) == 0) {
#        $currentregion = $region;
#      } else {
#        $currentregion .= "," . $region;
#      }
#
#      if ($hitEnd == 0) {
#        push(@retVal, $currentregion);
#        $currentregion = "";
#        $approxWindowRemaining = $approxWindow;
#      }
#    }
#  }
#
#  if (length($currentregion) > 0) {
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
  if (defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $nbJobs = LoadConfig::getParam($rH_cfg, 'mpileup', 'approxNbJobs');
  my $rA_regions = generateApproximateWindows($nbJobs, $rAoH_seqDictionary);

  my $bcfDir = "variants/rawBCF/";
  my $outputDir = "variants/";

  my $rO_job = SAMtools::mergeFilterBCF($rH_cfg, "allSamples", $bcfDir, $outputDir, $rA_regions);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeFilterBCF", undef, "MERGEBCF", $jobDependency, "allSamples", $rO_job);
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
  if (defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $inputVCF = "variants/allSamples.merged.flt.vcf";
  my $outputVCF = "variants/allSamples.merged.flt.NFiltered.vcf";

  my $rO_job = Tools::filterNStretches($rH_cfg, "allSamples", $inputVCF, $outputVCF);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "filterNStretches", undef, "FILTERN", $jobDependency, "allSamples", $rO_job);
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
  if (defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $inputVCF = "variants/allSamples.merged.flt.NFiltered.vcf";
  my $outputVCF = "variants/allSamples.merged.flt.mil.vcf";
  my $rO_job = VCFtools::annotateMappability($rH_cfg, $inputVCF, $outputVCF);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "flagMappability", undef, "MAPPABILITY", $jobDependency, "allSamples", $rO_job);
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
  if (defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $inputVCF = "variants/allSamples.merged.flt.mil.vcf";
  my $vcfOutput = "variants/allSamples.merged.flt.mil.snpId.vcf";

  my $rO_job = SnpEff::annotateDbSnp($rH_cfg, $inputVCF, $vcfOutput);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "snpIDAnnotation", undef, "SNPID", $jobDependency, "allSamples", $rO_job);
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
  if (defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }

  my $inputVCF = "variants/allSamples.merged.flt.mil.snpId.vcf";
  my $vcfOutput = "variants/allSamples.merged.flt.mil.snpId.snpeff.vcf";

  my $rO_job = SnpEff::computeEffects($rH_cfg, $inputVCF, $vcfOutput, 1);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "snpEffect", undef, "SNPEFF", $jobDependency, "allSamples", $rO_job);
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
  if (defined($globalDep{$parentStep}->{'experiment'})) {
    $jobDependency = $globalDep{$parentStep}->{'experiment'};
  }


  my $inputVCF = "variants/allSamples.merged.flt.mil.snpId.snpeff.vcf";
  my $vcfOutput = "variants/allSamples.merged.flt.mil.snpId.snpeff.dbnsfp.vcf";

  my $rO_job = SnpEff::annotateDbNSFP($rH_cfg, $inputVCF, $vcfOutput);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "dbNSFPAnnotation", undef, "DBNSFP", $jobDependency, "allSamples", $rO_job);
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

  if (defined($globalDep{$parentStep}->{'experiment'})) {
      $vcfDependency .= $globalDep{$parentStep}->{'experiment'};
  }


  my $inputVCF = "variants/allSamples.merged.flt.mil.snpId.vcf";
  my $outputFile = "metrics/allSamples.merged.flt.mil.snpId.snpeff.vcf.part.changeRate.tsv";
  my $listFiles = "variants/allSamples.merged.flt.mil.snpId.snpeff.vcf.statsFile.txt";

  my $rO_job_changeRate = Metrics::svnStatsChangeRate($rH_cfg, $inputVCF, $outputFile, $listFiles);
  if (!$rO_job_changeRate->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "metricsSNV", undef, "CHANGERATE", $vcfDependency, "allSamples", $rO_job_changeRate);
  }

  my $outputBaseName = "metrics/allSamples.SNV";
  my $rO_job_Graph = Metrics::svnStatsGetGraph($rH_cfg, $listFiles, $outputBaseName);
  if (!$rO_job_Graph->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "metricsSNV", "2", "SNV_GRAPH", $rO_job_changeRate->getCommandJobId(0), "allSamples", $rO_job_Graph);
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
    if (defined($globalDep{$stepName}->{'experiment'})) {
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $globalDep{$stepName}->{'experiment'};
    }
  }
  if (length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }
  $reportDependency = $jobDependencies;


  my $rO_job = GqSeqUtils::clientReport($rH_cfg, $configFile, $workDirectory, "DNAseq");

  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "deliverable", "REPORT", "DNAREPORT", $reportDependency, "allSamples", $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

1;
