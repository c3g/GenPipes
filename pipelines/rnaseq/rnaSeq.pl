#!/usr/bin/env perl

=head1 NAME

I<rnaSeq>

=head1 SYNOPSIS

  perl rnaSeq.pl -c rnaSeq.abacus.ini -s 1 -e 14 -n project.nanuq.csv -d design.txt -w  `pwd`  > toRun.sh

  will generate a bash script for steps 1 to 14. This script can then be executed:

  sh toRun.sh

  Options

  -c (rnaSeq.abacus.ini) the standard configuration file for the pipeline. Templates for some cluster systems like Abacus or Guillimin may already be available at pipelines/rnaseq/
  -s The start step
  -e The end step
  -n (project.nanuq.csv)  the NANUQ project read set sheet, prepared as described above.
  -d (design.txt) the design file. A tab separated value file that specifies the experimental design information of the project. The first column lists the sample names, which should match elements the column Name in the read set sheet. Subsequent columns specify all the pairwise comparisons which should be undertaken: values should be either "2" (nominator), "1" (denominator) or "0" (exclude from comparison).
  -w The project's working directory. All job outputs will be sent to this directory.


=head1 DESCRIPTION

B<rnaSeq> Is the main RNAseq pipeline.

=head1 AUTHOR

B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing

B<LoadConfig> Parse configuration file

B<Picard> Multiple tools to manage bam files (merge, sort, etc)

B<SampleSheet> Parse sample sheet file

B<SAMtools> alignment files (sam / bam) tools

B<SequenceDictionaryParser> Parse sequence dictionnary

B<SubmitToCluster> Create the submit command, control for dependencies

B<TophatBowtie> Alignment of RNA-Seq reads to a reference genome

B<Trimmomatic> Trim and filter raw read files

B<Metrics> Read, alignment and multiple metrics library

B<Cufflinks> Transcript assembly, differential expression, and differential regulation for RNA-Seq

B<Wiggle> Tools to generate wiggle tracks

B<HtseqCount> htseq is a library to generate basic statistics on aligned read count

B<DiffExpression> is a library to launch differential expression analysis (edgeR, deSeq, goseq)

B<GqSeqUtils>  is a library to access/launch functions from the gqSeqUtils R package

B<Version> Tracks app version


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
use Cwd qw/ abs_path /;
use File::Basename;
use Parse::Range qw(parse_range);

use LoadConfig;
use Picard;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use SubmitToCluster;
use TophatBowtie;
use Trimmomatic;
use Metrics;
use Cufflinks;
use Wiggle;
use HtseqCount;
use DiffExpression;
use GqSeqUtils;
use Version;

#--------------------


# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'samToFastq' , 'stepLoop' => 'sample' , 'output' => 'raw_reads'});
push(@steps, {'name' => 'trimming' , 'stepLoop' => 'sample' , 'output' => 'reads'});
push(@steps, {'name' => 'trimMetrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
push(@steps, {'name' => 'aligning' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
push(@steps, {'name' => 'merging' , 'stepLoop' => 'sample' , 'output' => 'alignment'});
push(@steps, {'name' => 'alignMetrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
#push(@steps, {'name' => 'mutation' , 'stepLoop' => 'sample' , 'output' => 'mpileup'});
push(@steps, {'name' => 'wiggle' , 'stepLoop' => 'sample' , 'output' => 'tracks'});
push(@steps, {'name' => 'rawCounts' , 'stepLoop' => 'sample' , 'output' => 'raw_counts'});
push(@steps, {'name' => 'rawCountsMetrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
push(@steps, {'name' => 'fpkm' , 'stepLoop' => 'sample' , 'output' => 'fpkm'});
push(@steps, {'name' => 'exploratory' , 'stepLoop' => 'group' , 'output' => 'exploratory'});
push(@steps, {'name' => 'cuffdiff' , 'stepLoop' => 'group' , 'output' => 'DGE'});
#push(@steps, {'name' => 'dgeMetrics' , 'stepLoop' => 'group' , 'output' => 'metrics'});
push(@steps, {'name' => 'dge' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'goseq' , 'stepLoop' => 'group' , 'output' => 'DGE'});
push(@steps, {'name' => 'deliverable' , 'stepLoop' => 'group' , 'output' => 'deliverable'});


#--------------------
# PODS
#--------------------
## Here starts the pipeline steps documentation, please change it accordingly any time you add/remove/modify a step

=head1 RNASEQ PIPELINE STEPS

The standard differential expression analysis for RNAseq data performs the following steps:

B<trimming> :  Raw reads quality trimming and removing of Illumina adapters is performed using trimmomatic.

B<trimMetrics> : Generates the trimming statistics file

B<aligning> : The filtered reads are aligned to a reference genome. The alignment is done per lane of sequencing using the combination of tophat/bowtie software. It generates a Binary Alignment Map file (.bam).

B<merging> : Bam files per sample are merged in one file. Merge is done using the Picard software. The resulting alignment file is reordered (karyotypic) and resulting reads per sample are marked as duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score) will be marked as a duplicate in the .bam file. Marking duplicates and reorder are executed using the Picard software.

B<wiggle> : generate wiggle tracks suitable for multiple browsers.

B<rawCounts> : Counting reads in features using htseq-count.

B<rawCountsMetrics> : Create rawcount matrix, zip the wiggle tracks and create the saturation plots based on standardized read counts.

B<fpkm> : Compute FPKM measures for de novo and known transcripts using cufflinks.

B<exploratory> : exploratoryAnalysis using the gqSeqUtils R package

B<cuffdiff> : The transcript quantification engine of Cufflinks (Cuffdiff) is used to calculate transcript expression levels in more than one condition and test them for signficant differences.

B<dge> : Differential gene expression analysis using DESEQ and EDGER. Merge the results of the analysis in a single csv file.

B<goseq>:  Gene Ontology analysis for RNA-seq using the bioconductor's R package goseq. Generates Go annotations for cuffdiff known transcripts differential expression and differential gene expression analysis.

B<deliverable> : Generates the standard report. A summary html report contains the description of the sequencing experiment as well as a detailed presentation of the pipeline steps and results. Various Quality Control (QC) summary statistics are included in the report and additional QC analysis is accessible for download directly through the report. The report includes also the main references of the software and methods used during the analysis, together with the full list of parameters passed to the pipeline main script.

=cut end of documentation


my %globalDep;
for my $stepName (@steps) {
  $globalDep{$stepName->{'name'} } = {};
}


# Global scope variables
my $designFilePath;
my $configFile;
my $workDir;
my $readSetSheet;


&main();

sub printUsage {
  print "Version: " . $Version::version . "\n";
  print "\nUsage: perl " . $0 . " \n";
  print "\t-c  config file\n";
  print "\t-s  step range e.g. '1,3', '2-5', '1,4-7,10'\n";
  print "\t-n  nanuq sample sheet\n";
  print "\t-d  design file\n";
  print "\t-w  work directory\n";
  print "\n";
  print "Steps:\n";
  for (my $idx = 0; $idx < @steps; $idx++) {
    print "" . ($idx + 1) . '- ' . $steps[$idx]->{'name'} . "\n";
  }
  print "\n";
}

sub main {
	my %opts;
	getopts('c:s:e:n:d:w:', \%opts);
	
	if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'n'}) || !defined($opts{'w'} ) ) {
		printUsage();
		exit(1);
	}
	
	my %jobIdVarPrefix;
	my %cfg = LoadConfig->readConfigFile($opts{'c'});
	my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($opts{'n'});
	my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);
  
	$designFilePath = $opts{'d'};
  my $rHoAoA_designGroup;
  if (defined($designFilePath)) {
    $designFilePath = abs_path($designFilePath);
	  ##get design groups
	  $rHoAoA_designGroup = Cufflinks::getDesign(\%cfg,$designFilePath);
  }
	$workDir = abs_path($opts{'w'});
	$configFile =  abs_path($opts{'c'});
	$readSetSheet =  abs_path($opts{'n'}); 

	
	#generate sample jobIdprefix
	my $cpt = 1;
	
	for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
		my $cpt2=1;
		$jobIdVarPrefix{$sampleName} = $cpt;
		my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
		for my $rH_laneInfo (@$rAoH_sampleLanes) {
			$jobIdVarPrefix{$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} = $cpt.'_'.$cpt2;
			$cpt2++;
		}
		$cpt++;
	}
	#generate design jobIdprefix
	if(defined($rHoAoA_designGroup)) {
  	for my $designName (keys %{$rHoAoA_designGroup}) {
	  	$jobIdVarPrefix{$designName} = $cpt;
		  $cpt++;
  	}
  }
	
  # List user-defined step index range.
  # Shift 1st position to 0 instead of 1
  my @stepRange = map($_ - 1, parse_range($opts{'s'}));

	SubmitToCluster::initPipeline($workDir);

	for my $current (@stepRange) {
		my $fname = $steps[$current]->{'name'};
		my $loopType = $steps[$current]->{'stepLoop'};
		my $subref = \&$fname;
		if ($loopType eq 'sample') {
			for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
			my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
			# Tests for the first step in the list. Used for dependencies.
			my $jobIdVar = &$subref($current != $stepRange[0], \%cfg, $sampleName, $rAoH_sampleLanes, $rAoH_seqDictionary, \%jobIdVarPrefix);
			if (defined($jobIdVar)) {
				$globalDep{$fname}->{$sampleName} = $jobIdVar;
			}
			}
		}
		else {
			# Tests for the first step in the list. Used for dependencies.
			my $jobIdVar = &$subref($current != $stepRange[0], \%cfg, $rHoAoH_sampleInfo, $rHoAoA_designGroup, $rAoH_seqDictionary, \%jobIdVarPrefix);
			if (defined($jobIdVar)) {
				$globalDep{$fname}->{$fname} = $jobIdVar;
			}
		}
	}
}

sub samToFastq {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $samToFastqJobIdVarNameSample = undef;
	my $libraryType = LoadConfig::getParam($rH_cfg, 'default', 'libraryType');
	for my $rH_laneInfo (@$rAoH_sampleLanes) {

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
        die "Error in rnaSeqDeNovoAssembly::samToFastq: unknown run type (can be 'SINGLE_END' or 'PAIRED_END' only): " . $rH_laneInfo->{'runType'};
      }
    }

		if(!$rO_job->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "samToFastq", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'SAMTOFASTQ' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $rO_job);
      if(!defined($samToFastqJobIdVarNameSample)) {
			  $samToFastqJobIdVarNameSample = $rO_job->getCommandJobId(0);
      } else {
        $samToFastqJobIdVarNameSample .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
      }
		}
	}
	return $samToFastqJobIdVarNameSample;	
}

sub trimming {
	my $depends = shift;
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rAoH_sampleLanes  = shift;
	my $rAoH_seqDictionary = shift;
	my $rH_jobIdPrefixe = shift;

	my $trimJobIdVarNameSample = undef;
	my $libraryType = LoadConfig::getParam($rH_cfg, 'default', 'libraryType');
	my $jobDependency = undef;
	if($depends > 0) {
	  $jobDependency = $globalDep{'samToFastq'}{$sampleName};
	}
	for my $rH_laneInfo (@$rAoH_sampleLanes) {
		#print "mkdir -p metrics/$sampleName/output_jobs reads/$sampleName/output_jobs\n";
		##get raw read count
# 		my $inputFile = LoadConfig::getParam($rH_cfg, 'default', 'rawReadDir', 1, 'dirpath') .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' .$rH_laneInfo->{'read1File'};
# 		my $outputFile= 'metrics/' .$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.raw.csv' ;
# 		my $command = Metrics::readStats($rH_cfg,$inputFile,$outputFile,'fastq',$libraryType);
# 		my $rawReadStatJobID = undef;
# 		if(defined($command) && length($command) > 0) {
# 			$rawReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'raw', 'RAWREADSTAT' .$rH_jobIdPrefixe ->{$sampleName .'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $command);
# 			$rawReadStatJobID = '$'.$rawReadStatJobID;
# 		}
# 		
		## trimming - TO DO should be modified to the new rawread location (cf. David modif) and portability
		my $minQuality  = $rH_cfg->{'trim.minQuality'};
		my $minLength   = $rH_cfg->{'trim.minLength'};
		my $laneDir = 'reads/' .$sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
		print "mkdir -p $laneDir\n";
		my $outputFastqPair1Name;
		if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
			$outputFastqPair1Name = $laneDir . $sampleName.'.'.$rH_laneInfo->{'libraryBarcode'}.'.t'.$minQuality.'l'.$minLength.'.single.fastq.gz';
		}
		elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
			$outputFastqPair1Name = $laneDir . $sampleName.'.'.$rH_laneInfo->{'libraryBarcode'}.'.t'.$minQuality.'l'.$minLength.'.pair1.fastq.gz';
		}
		else {
			die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
		}

		my $rO_job = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo, $laneDir);
		if(!$rO_job->isUp2Date()) {
			SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM' .$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , $jobDependency, $sampleName, $rO_job);
      if(!defined($trimJobIdVarNameSample)) {
			  $trimJobIdVarNameSample = $rO_job->getCommandJobId(0);
      } else {
        $trimJobIdVarNameSample .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
      }
    }
  }
  return $trimJobIdVarNameSample;
}


sub trimMetrics {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $libraryType = LoadConfig::getParam($rH_cfg, 'default', 'libraryType');
  my $trimmingDependency = undef;
  if ($depends > 0) {
    $trimmingDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), values(%{$globalDep{'trimming'}}));
  }
  print "mkdir -p metrics/\n";
  my $folder = 'reads';
  my $pattern = 'trim.stats.csv';
  my $ouputFile = 'metrics/trimming.stats';
  my $rO_job = Metrics::mergeTrimmomaticStats($rH_cfg, $libraryType, $pattern, $folder, $ouputFile);
  my $metricsJobId = undef;
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "trimMetrics", undef, 'TRIMMETRICS', $trimmingDependency, undef, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}

sub aligning {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $alignJobIdVarNameSample = undef;
  my $jobDependency = undef;
  if ($depends > 0) {
    $jobDependency = $globalDep{'trimming'}{$sampleName};
  }

  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $pair1;
    my $pair2;
    my $single;
    #align lanes
    my $outputDirPath = 'alignment/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    print "mkdir -p $outputDirPath \n";
    my $rO_job;
    if ($rH_laneInfo->{'runType'} eq "SINGLE_END") {
      $single = 'reads/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int') . 'l' . LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int') . '.single.fastq.gz';
      $rO_job = TophatBowtie::align($rH_cfg, $sampleName, $rH_laneInfo, $single, undef);
    } elsif ($rH_laneInfo->{'runType'} eq "PAIRED_END") {
      $pair1 = 'reads/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int') . 'l' . LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int') . '.pair1.fastq.gz';
      $pair2 = 'reads/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName . '.' . $rH_laneInfo->{'libraryBarcode'} . '.t' . LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int') . 'l' . LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int') . '.pair2.fastq.gz';
      $rO_job = TophatBowtie::align($rH_cfg, $sampleName, $rH_laneInfo, $pair1, $pair2);
    }

    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "align", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'ALIGN' . $rH_jobIdPrefix ->{$sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} . 'ALIGN', $jobDependency, $sampleName, $rO_job);
      if (!defined($alignJobIdVarNameSample)) {
        $alignJobIdVarNameSample = $rO_job->getCommandJobId(0);
      } else {
        $alignJobIdVarNameSample .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
      }
    }
    ##generate aligment stats
#     my $inputFile = 'alignment/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . 'accepted_hits.bam';
#     my $outputFile= 'metrics/' . $sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.aligned.csv';
#     $command = Metrics::readStats($rH_cfg, $inputFile, $outputFile, 'bam');
#     my $alignedReadStatJobID = undef;
#     if (defined($command) && length($command) > 0) {
#       $alignedReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'aligned', 'ALIGNEDREADSTAT' . $rH_jobIdPrefix->{$sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}}, $alignJobIdVarNameLane, $sampleName, $command);
#       $alignedReadStatJobID = '$' . $alignedReadStatJobID;
#     }
#     ##merge read stats
#     my $rawFile = 'metrics/' . $sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.raw.csv';
#     my $filterFile = 'metrics/' . $sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.filtered.csv';
#     my $alignFile = $outputFile;
#     my $sampleNameFull = $sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
#     $outputFile= 'metrics/' . $sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '.readstats.csv';
#     $command = Metrics::mergeIndvidualReadStats($rH_cfg, $sampleNameFull, $rawFile, $filterFile, $alignFile, $outputFile);
#     my $mergeReadStatJobID = undef;
#     if (defined($command) && length($command) > 0) {
#       $mergeReadStatJobID = SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'merged', 'MERGEREADSTAT' . $rH_jobIdPrefix->{$sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}}, $alignedReadStatJobID, $sampleName, $command);
#       $mergeReadStatJobID = '$' . $mergeReadStatJobID;
#       $alignJobIdVarNameSample .= $mergeReadStatJobID . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
#     }
  }
  return $alignJobIdVarNameSample;
}

sub merging {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $jobDependency = undef;
  if ($depends > 0) {
    $jobDependency = $globalDep{'aligning'}{$sampleName};
  }

  ##Merging
  my $inputBAM;
  my $outputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.bam';
  my @alignFiles;
  my $merge = 0; # JT : Flag to see if we skip merging step or not.

  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $laneDir = "alignment/" . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    $inputBAM = $laneDir . 'accepted_hits.bam';
    $merge++; # JT: increment flag
    push(@alignFiles, $inputBAM);
  }

  my $mergeJobId;
  if ($merge > 1) { # JT: If flag is higher than 1 (so more than one lane / sample), perform merge.
    my $rO_mergeJob = Picard::mergeFiles($rH_cfg, $sampleName, \@alignFiles, $outputBAM);
    if (!$rO_mergeJob->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "mergeFiles", undef, 'MERGELANES' . $rH_jobIdPrefix ->{$sampleName}, $jobDependency, $sampleName, $rO_mergeJob);
      $mergeJobId = $rO_mergeJob->getCommandJobId(0);
    }
  } else {
    $mergeJobId = $jobDependency; # JT : update job dependency as well.
  }

  ## reorder
  if ($merge > 1) { # JT : if merge, update input file accordingly.
    $inputBAM = $outputBAM; #Update file name if it has been merged.
  }
  $outputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.karyotypic.bam';
  my $rO_reorderJob = Picard::reorderSam($rH_cfg, $sampleName, $inputBAM, $outputBAM);
  if (!$rO_reorderJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "reorderSam", undef, 'REORDER' . $rH_jobIdPrefix->{$sampleName} . 'REORDER', $mergeJobId, $sampleName, $rO_reorderJob);
  }

  ## mark duplicates
  $inputBAM = $outputBAM;
  $outputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam';
  my $duplicatesMetricsFile = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.metrics';
  my $rO_markDupJob = Picard::markDup($rH_cfg, $sampleName, $inputBAM, $outputBAM, $duplicatesMetricsFile);
  if (!$rO_markDupJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'MARKDUP' . $rH_jobIdPrefix->{$sampleName}, $rO_reorderJob->getCommandJobId(0), $sampleName, $rO_markDupJob);
  }
  return $rO_markDupJob->getCommandJobId(0);
}



sub alignMetrics {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $mergingDependency = undef;
  if ($depends > 0) {
    $mergingDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), values(%{$globalDep{'merging'}}));
  }
  ## RNAseQC metrics
  mkdir $workDir . '/alignment';
  my $sampleList = $workDir . '/alignment/rnaseqc.samples.txt';
  open(RNASAMPLE, ">$sampleList") or die ("Unable to open $sampleList for writing");
  print RNASAMPLE "Sample\tBamFile\tNote\n";
  my $projectName = LoadConfig::getParam($rH_cfg, 'metricsRNA', 'projectName');
  for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
    print RNASAMPLE "$sampleName\talignment/$sampleName/$sampleName.merged.mdup.bam\t$projectName\n";
  }
  close(RNASAMPLE);
  print "mkdir -p metrics/rnaseqRep/\n";
  $sampleList = 'alignment/rnaseqc.samples.txt';
  my $outputFolder = 'metrics/rnaseqRep';
  my $libraryType = LoadConfig::getParam($rH_cfg, 'default', 'libraryType');
  my $rO_job = Metrics::rnaQc($rH_cfg, $sampleList, $outputFolder, $libraryType);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "rnaQc", LoadConfig::getParam($rH_cfg, 'rnaQc', 'projectName'), 'METRICSRNA', $mergingDependency, undef, $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}


# sub mutation {
#   my $depends = shift;
#   my $rH_cfg = shift;
#   my $sampleName = shift;
#   my $rAoH_sampleLanes = shift;
#   my $rAoH_seqDictionary = shift;
# ALIGNEDREADSTAT_JOB_ID
#   my $jobDependency = undef;
#   if ($depends > 0) {
#     $jobDependency = $globalDep{'merging'}{$sampleName};
#   }
#
#   my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.bam";
#
#
# }


sub wiggle {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $jobDependency = undef;
  if ($depends > 0) {
    $jobDependency = $globalDep{'merging'}{$sampleName};
  }

  my $inputBAM = "alignment/" . $sampleName . "/" . $sampleName . ".merged.mdup.bam";
  #testing for strand-specificity

  my $strandSPecificityInfo = LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
  my @strandJobId;
  my @outputBAM;
  my @outputBedGraph;
  my @outputWiggle;
  my @prefixJobName;
  print "mkdir -p tracks/$sampleName/ tracks/bigWig/\n";
  if ($strandSPecificityInfo ne "fr-unstranded") {
  ## strand specific
    @outputBAM = ('alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.forward.bam', 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.reverse.bam');
    my $rO_job = Wiggle::strandBam($rH_cfg, $sampleName, $inputBAM, \@outputBAM);
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", 'FORWARD', 'FSTRANDSPEC' . $rH_jobIdPrefix->{$sampleName}, $jobDependency, $sampleName, $rO_job, 0);
      push(@strandJobId, $rO_job->getCommandJobId(0));
      SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", 'REVERSE', 'RSTRANDSPEC' . $rH_jobIdPrefix->{$sampleName}, $jobDependency, $sampleName, $rO_job, 1);
      push(@strandJobId, $rO_job->getCommandJobId(1));
    }
    @outputBedGraph = ('tracks/' . $sampleName . '/' . $sampleName . '.forward.bedGraph', 'tracks/' . $sampleName . '/' . $sampleName . '.reverse.bedGraph');
    @outputWiggle = ('tracks/bigWig/' . $sampleName . '.forward.bw', 'tracks/bigWig/' . $sampleName . '.reverse.bw');
    @prefixJobName = ('FORWARD2', 'REVERSE2');
  } else {
    push(@outputBAM, $inputBAM);
    push(@strandJobId, $jobDependency);
    push(@outputBedGraph, 'tracks/' . $sampleName . '/' . $sampleName . '.bedGraph');
    push(@outputWiggle, 'tracks/bigWig/' . $sampleName . '.bw');
    push(@prefixJobName, undef);
  }
  my $wiggleJobId = undef;
  for (my $i = 0; $i <@outputBAM; $i++) {
    my $rO_job = Wiggle::graph($rH_cfg, $sampleName, $inputBAM, $outputBedGraph[$i], $outputWiggle[$i]);
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", $prefixJobName[$i], 'WIGGLE' . $rH_jobIdPrefix->{$sampleName}, $strandJobId[$i], $sampleName, $rO_job);
      if (!defined($wiggleJobId)) {
        $wiggleJobId = $rO_job->getCommandJobId(0);
      } else {
        $wiggleJobId .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
      }
    }
  }
  return $wiggleJobId;
}


sub rawCounts {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $jobDependency = undef;
  if ($depends > 0) {
    $jobDependency = $globalDep{'merging'}{$sampleName};
  }
  print "mkdir -p raw_counts\n";
  my $inputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam';
  my $sortedBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.queryNameSorted.bam';
  my $inputGtf = LoadConfig::getParam($rH_cfg, 'htseq', 'referenceGtf', 1, 'filepath');
  my $outputCount = 'raw_counts/' . $sampleName . '.readcounts.csv';
  my $sortOrder = 'queryname';
  my $strandInfo;
  my $strandSPecificityInfo = LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
  if ($strandSPecificityInfo ne "fr-unstranded") {
    $strandInfo= 'reverse';
  } else {
    $strandInfo= 'no';
  }
  ## query sort the bam
  my $rO_sortJob = Picard::sortSam($rH_cfg, $sampleName, $inputBAM, $sortedBAM, $sortOrder);
  if (!$rO_sortJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "sortSam", undef, 'QNSORT' . $rH_jobIdPrefix->{$sampleName}, $jobDependency, $sampleName, $rO_sortJob);
  }
  ## count reads
  my $rO_countJob = HtseqCount::readCountPortable($rH_cfg, $sortedBAM, $inputGtf, $outputCount, $strandInfo);
  if (!$rO_countJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "htseq", undef, 'RAWCOUNT' . $rH_jobIdPrefix->{$sampleName}, $rO_sortJob->getCommandJobId(0), $sampleName, $rO_countJob);
  }
  return $rO_countJob->getCommandJobId(0);
}


sub rawCountsMetrics {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $metricsJobId = undef;
  my $countDependency = undef;
  my $wiggleDependency = undef;
  if ($depends > 0) {
    $countDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), values(%{$globalDep{'rawCounts'}}));
    $wiggleDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), values(%{$globalDep{'wiggle'}}));
  }
  #create rawcount matrix
  print "mkdir -p DGE\n";
  my $readCountDir = 'raw_counts';
  my $readcountExtension = '.readcounts.csv';
  my $outputDir = 'DGE';
  my $outputMatrix = 'rawCountMatrix.csv';
  my $rO_matrixJob = HtseqCount::refGtf2matrix($rH_cfg, LoadConfig::getParam($rH_cfg, 'htseq', 'referenceGtf', 1, 'filepath'), $readCountDir, $readcountExtension, $outputDir, $outputMatrix);
  if (!$rO_matrixJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", 'matrix', 'MATRIX', $countDependency, undef, $rO_matrixJob);
    $metricsJobId = $rO_matrixJob->getCommandJobId(0);
  }

  ### to do outside of the wiggle function on time only
  print "mkdir -p tracks/bigWig\n";
  my $wigFolder = 'tracks/bigWig/';
  my $wigArchive = 'tracks.zip';
  my $rO_job = Wiggle::zipWig($rH_cfg, $wigFolder, $wigArchive);
  my $wigZipJobId;
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "metrics", undef, 'WIGZIP', $wiggleDependency, undef, $rO_job);
    if (!defined($metricsJobId)) {
      $metricsJobId = $rO_job->getCommandJobId(0);
    } else {
      $metricsJobId .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
    }
  }
  ##RPKM and Saturation
  print "mkdir -p metrics/saturation\n";;
  my $countFile = 'DGE/rawCountMatrix.csv';
  my $geneSizeFile = LoadConfig::getParam($rH_cfg, 'saturation', 'geneSizeFile', 1, 'filepath');
  my $rpkmDir = 'raw_counts';
  my $saturationDir = 'metrics/saturation';

  $rO_job = Metrics::saturation($rH_cfg, $countFile, $geneSizeFile, $rpkmDir, $saturationDir);
  my $saturationJobId = undef;
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "saturation", undef, 'RPKM', $rO_matrixJob->getCommandJobId(0), undef, $rO_job);
    if (!defined($metricsJobId)) {
      $metricsJobId = $rO_job->getCommandJobId(0);
    } else {
      $metricsJobId .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
    }
  }
  return $metricsJobId;
}

sub fpkm {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $jobDependency = undef;
  if ($depends > 0) {
    $jobDependency = $globalDep{'merging'}{$sampleName};
  }
  print "mkdir -p fpkm/known/$sampleName fpkm/denovo/$sampleName\n";
  my $inputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.mdup.bam';
  my $outputKnown = 'fpkm/known/' . $sampleName;
  my $outputDeNovo = 'fpkm/denovo/' . $sampleName;
  my $gtfOption = '-G ' . LoadConfig::getParam($rH_cfg, 'fpkm', 'referenceGtf', 1, 'filepath');

  ## known FPKM
  my $fpkmJobId = undef;
  my $rO_fpkmKnownJob = Cufflinks::fpkm($rH_cfg, $inputBAM, $outputKnown, $gtfOption);
  if (!$rO_fpkmKnownJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "fpkm", "KNOWN", 'FPKMK' . $rH_jobIdPrefix->{$sampleName} . 'FPKM', $jobDependency, $sampleName, $rO_fpkmKnownJob);
    $fpkmJobId = $rO_fpkmKnownJob->getCommandJobId(0);
  }
  ## denovo FPKM
  my $rO_fpkmDeNovoJob = Cufflinks::fpkm($rH_cfg, $inputBAM, $outputDeNovo, undef);
  if (!$rO_fpkmDeNovoJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "fpkm", "DENOVO", 'FPKMD' . $rH_jobIdPrefix->{$sampleName}, $jobDependency, $sampleName, $rO_fpkmDeNovoJob);
    if (!defined($fpkmJobId)) {
      $fpkmJobId = $rO_fpkmDeNovoJob->getCommandJobId(0);
    } else {
      $fpkmJobId .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_fpkmDeNovoJob->getCommandJobId(0);
    }
  }
  return $fpkmJobId;
}


sub cuffdiff {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $jobDependency = undef;
  if ($depends > 0 and values(%{$globalDep{'fpkm'}}) > 0) {
    $jobDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), values(%{$globalDep{'fpkm'}}));
  }

  if (!defined($rHoAoA_designGroup)) {
    die ("A design is needed for the cuffdiff step\n");
  }

  print "mkdir -p cuffdiff/known cuffdiff/denovo\n";
  ## create the list of deNovo gtf to merge
  mkdir $workDir . '/fpkm';
  mkdir $workDir . '/fpkm/denovo/';
  my $mergeListFile = $workDir . '/fpkm/denovo/gtfMerge.list';
  my @compareList;
  open(MERGEF, ">$mergeListFile") or die ("Unable to open $mergeListFile for wrtting");
  ##iterate over sample
  for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
    my $gtfFile = 'fpkm/denovo/' . $sampleName . '/transcripts.gtf';
    push(@compareList, 'fpkm/denovo/' . $sampleName . '/transcripts.gtf');
    print MERGEF $gtfFile;
    print MERGEF "\n";
  }
  close($mergeListFile);
  ##merge denovo transcript in one gtf file
  my $outputPathDeNovo = 'fpkm/denovo/allSample';
  my $rO_job = Cufflinks::cuffcompare($rH_cfg, \@compareList, $outputPathDeNovo, $mergeListFile);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "cuffcompare", "MERGE", 'GTFCOMPARE', $jobDependency, undef, $rO_job);
  }
#  my $gtfDnMerged = 'fpkm/denovo/merged.gtf';
#  my $gtfDnFormatMerged = 'fpkm/denovo/formated.merged.gtf';
#  $command = Cufflinks::mergeGtfFormat($rH_cfg, $gtfDnMerged, $gtfDnFormatMerged);
#  my $formatJobId;
#  if (defined($command) && length($command) > 0) {
#    $formatJobId= SubmitToCluster::printSubmitCmd($rH_cfg, "default", "FORMAT", 'GTFFORMAT', $cuffmergeJobId, undef, $command);
#    $formatJobId= '$' . $formatJobId;
#    $cuffJobID = $formatJobId . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
#  }
  ##iterate over design
  my $cuffddiffJobId = undef;
  for my $design (keys %{$rHoAoA_designGroup}) {
    mkdir $workDir;
    mkdir $workDir . '/cuffdiff';
    my $numberGroups = @{$rHoAoA_designGroup->{$design}};
    ##iterate over group
    my @groupInuptFiles;
    for (my $i = 0; $i < $numberGroups; $i++) {
      ##iterate over samples in the design
      my $numberSample = @{$rHoAoA_designGroup->{$design}->[$i]};
      my $bamfile = ' ';
      for (my $j = 0; $j < $numberSample; $j++) {
        $bamfile .= 'alignment/' . $rHoAoA_designGroup->{$design}->[$i]->[$j] . '/' . $rHoAoA_designGroup->{$design}->[$i]->[$j] . '.merged.mdup.bam' .',';
      }
      $bamfile = substr $bamfile, 0, -1;
      push(@groupInuptFiles, $bamfile);
    }

    my $outputPathKnown = 'cuffdiff/known/' . $design;
    ##cuffdiff known
    $rO_job = Cufflinks::cuffdiff($rH_cfg, \@groupInuptFiles, $outputPathKnown, LoadConfig::getParam($rH_cfg, 'cuffdiff', 'referenceGtf', 1, 'filepath'));
    if (!$rO_job->isUp2Date()) {
      my $diffJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "cuffdiff", "KNOWN", 'CUFFDIFFK' . $rH_jobIdPrefix->{$design}, $jobDependency, $design, $rO_job);
      if (!defined($cuffddiffJobId)) {
        $cuffddiffJobId = $rO_job->getCommandJobId(0);
      } else {
        $cuffddiffJobId .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
      }
    }
  }
  $rO_job = Cufflinks::mergeCuffdiffRes($rH_cfg, $designFilePath, 'cuffdiff', 'fpkm');
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "default", "MERGE_RES", 'CUFF_MERGE_RES', $cuffddiffJobId, undef, $rO_job);
  }
  return $rO_job->getCommandJobId(0);
}


sub dge {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $jobDependency = undef;
  if ($depends > 0) {
    $jobDependency = $globalDep{'rawCountsMetrics'}{'rawCountsMetrics'};
  }

  my $countMatrix = 'DGE/rawCountMatrix.csv';
  my $outputDir = 'DGE';

  if (!defined($rHoAoA_designGroup)) {
    die ("A design is needed for the dge step\n");
  }

  ## edgeR
  my $rO_job = DiffExpression::edgerPortable($rH_cfg, $designFilePath, $countMatrix, $outputDir);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "diffExpress", 'edger', 'EDGER', $jobDependency, undef, $rO_job);
  }

  ## DESeq
  my $rO_deseqJob = DiffExpression::deseq($rH_cfg, $designFilePath, $countMatrix, $outputDir);
  if (!$rO_deseqJob->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "diffExpress", 'deseq', 'DESEQ', $rO_job->getCommandJobId(0), undef, $rO_deseqJob);
  }

  return $rO_deseqJob->getCommandJobId(0);
}

sub goseq {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $cuffdiffDependency = undef;
  my $dgeDependency = undef;
  if ($depends > 0) {
    $cuffdiffDependency = $globalDep{'cuffdiff'}{'cuffdiff'};
    $dgeDependency = $globalDep{'dge'}{'dge'};
  }

  if (!defined($rHoAoA_designGroup)) {
    die ("A design is needed for the goseq step\n");
  }

  my $columnsCuff = LoadConfig::getParam($rH_cfg, 'diffExpress', 'cuffRescolumns');
  my $columnsDge = LoadConfig::getParam($rH_cfg, 'diffExpress', 'dgeRescolumns');
  my $goseqJobId = undef;
  for my $design (keys %{$rHoAoA_designGroup}) {
    ## goseq for cuffdiff known results
    my $resultFileCuff = 'cuffdiff/known/' . $design . '/isoform_exp.diff';
    my $outputFileCuff = 'cuffdiff/known/' . $design . '/gene_ontology_results.csv';
    my $rO_job = DiffExpression::goseq($rH_cfg, $resultFileCuff, $outputFileCuff, $columnsCuff);
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "goseq", "CUFFLINKS", 'GOCUFFDIFF' . $rH_jobIdPrefix->{$design}, $cuffdiffDependency, $design, $rO_job);
      if (!defined($goseqJobId)) {
        $goseqJobId = $rO_job->getCommandJobId(0);
      } else {
        $goseqJobId .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
      }
    }
    ## goseq for dge results
    my $resultFileDge = 'DGE/' . $design . '/dge_results.csv';
    my $outputFileDge = 'DGE/' . $design . '/gene_ontology_results.csv';
    $rO_job = DiffExpression::goseq($rH_cfg, $resultFileDge, $outputFileDge, $columnsDge);
    if (!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "goseq", "DESEQ", 'GODGE' . $rH_jobIdPrefix->{$design}, $dgeDependency, $design, $rO_job);
      if (!defined($goseqJobId)) {
        $goseqJobId = $rO_job->getCommandJobId(0);
      } else {
        $goseqJobId .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $rO_job->getCommandJobId(0);
      }
    }
  }
  return $goseqJobId;
}

sub exploratory {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $jobDependency = undef;
  if ($depends > 0 and values(%{$globalDep{'fpkm'}}) > 0) {
    $jobDependency = join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), values(%{$globalDep{'fpkm'}}));
    if ($depends > 0 and values(%{$globalDep{'rawCountsMetrics'}}) > 0) {
      $jobDependency .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $globalDep{'rawCountsMetrics'}{'rawCountsMetrics'};
    }
  }

  print "mkdir -p exploratory\n";

  # Call gqSeqUtils::exploratoryAnalysis()
  # sub exploratoryRnaAnalysis {
  #         my $rH_cfg        = shift;
  #         my $readSetSheet  = shift;
  #         my $workDir       = shift;
  #         my $configFile    = shift;
  #
  my $rO_job = GqSeqUtils::exploratoryRnaAnalysis($rH_cfg, $readSetSheet, $workDir, $configFile);
  if (!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "exploratory", 'exploratoryAnalysis', 'exploratoryAnalysis', $jobDependency, undef, $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

sub deliverable {
  my $depends = shift;
  my $rH_cfg = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rHoAoA_designGroup = shift;
  my $rAoH_seqDictionary = shift;
  my $rH_jobIdPrefix = shift;

  my $jobDependency;
  if ($depends > 0) {
    my $trimDependency = $globalDep{'trimMetrics'}{'trimMetrics'};
    if (defined($trimDependency) && length($trimDependency) > 0) {
      $jobDependency .= $trimDependency . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
    }
    my $alignDependency = $globalDep{'alignMetrics'}{'alignMetrics'};
    if (defined($alignDependency) && length($alignDependency) > 0) {
      $jobDependency .= $alignDependency . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
    }
    my $goDependency = $globalDep{'goseq'}{'goseq'};
    if (defined($goDependency) && length($goDependency) > 0) {
      $jobDependency .= $goDependency . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
    }
  }

  if (defined($jobDependency) && length($jobDependency) > 0) {
    $jobDependency = substr $jobDependency, 0, -1;
  }

  my $rO_job = GqSeqUtils::clientReport($rH_cfg, $configFile, $workDir, 'RNAseq');
  my $deliverableJobId = undef;
  if (!$rO_job->isUp2Date()) {
    print "mkdir -p deliverable\n";
    SubmitToCluster::printSubmitCmd($rH_cfg, "deliverable", 'REPORT', 'RNAREPORT', $jobDependency, undef, $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

1;
