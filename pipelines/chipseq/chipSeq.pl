#!/usr/bin/env perl

=head1 NAME

I<chipSeq>

=head1 SYNOPSIS

  perl chipSeq.pl -c chipSeq.abacus.ini -n project.nanuq.csv -d design.csv -w  currentDirectory -s 1 -e 11 > toRun.sh

  Options:

    -c (chipSeq.abacus.ini) the standard configuration file for the pipeline.
    -s The start step
    -e The end step
    -n (project.nanuq.csv) the NANUQ Project sample file
    -d (design.csv) the design file. A tab separated value file that specifies the experimental design information of the project.
    -w The project's working directory. All job outputs will be sent to this directory.
  

=head1 DESCRIPTION

B<chipSeq.pl> is the main ChIPseq pipeline.

=head1 AUTHOR

B<Johanna Sandoval> - I<johanna.sandoval@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing

B<LoadConfig> Parse configuration file

B<SampleSheet> Parse sample sheet file

B<SAMtools> Samtools functions

B<SequenceDictionaryParser> Parse Picard generated dictionnary

B<Picard> Picard fonctions on bam files (merge, sort, etc)

B<SubmitToCluster> Paste submit options to the bash command, controls resume jobs

B<Homer>  Homer fonctions (Chipseq QC, annotation, motifs)

B<MACS2>  MACS peak call, design file validation

B<BWA> BWA alignment fonctions

B<Trimmomatic>  Trimmomatic trimming / clipping functions

B<Metrics>   Multiple metrics functions (trim/align/annotation)

B<Version>  Tracks app version

B<Wiggle >  Optional functions to generate wiggle tracks (deprecated)

B<GqSeqUtils>  Deliverable and final report generation


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
use Parse::Range qw(parse_range);
use POSIX;

use LoadConfig;
use SampleSheet;
use SAMtools;
use SequenceDictionaryParser;
use Picard;
use SubmitToCluster;
use Homer;
use MACS2;
use BWA;
use Trimmomatic;
use Metrics;
use Version;
use Wiggle;
use GqSeqUtils;



#--------------------
# SUB
#--------------------

my @steps;
push(@steps, {'name' => 'bamToFastq', 'stepLoop' => 'sample', 'output' => 'reads', 'parentStep' => undef});
push(@steps, {'name' => 'trimming', 'stepLoop' => 'sample', 'output' => 'reads', 'parentStep' => 'bamToFastq'});
push(@steps, {'name' => 'aligning', 'stepLoop' => 'sample', 'output' => 'alignment', 'parentStep' => 'trimming'});
push(@steps, {'name' => 'metrics', 'stepLoop' => 'sample', 'output' => 'metrics', 'parentStep' => 'aligning'});
push(@steps, {'name' => 'qcTagDirectories', 'stepLoop' => 'sample', 'output' => 'tags', 'parentStep' => 'aligning'});
push(@steps, {'name' => 'qcPlots', 'stepLoop' => 'experiment', 'output' => 'graphs', 'parentStep' => ['qcTagDirectories']});
push(@steps, {'name' => 'wiggle', 'stepLoop' => 'sample', 'output' => 'tracks', 'parentStep' => 'qcTagDirectories'});
push(@steps, {'name' => 'peakCall', 'stepLoop' => 'design', 'output' => 'peak_call', 'parentStep' => 'aligning'});
push(@steps, {'name' => 'annotation', 'stepLoop' => 'design', 'output' => 'annotation', 'parentStep' => 'peakCall'});
push(@steps, {'name' => 'motif', 'stepLoop' => 'design', 'output' => 'annotation', 'parentStep' => 'peakCall'});
push(@steps, {'name' => 'annotationPlots', 'stepLoop' => 'experiment', 'output' => 'graphs', 'parentStep' => ['annotation']});
push(@steps, {'name' => 'deliverable', 'stepLoop' => 'experiment', 'output' => 'Deliverable', 'parentStep' => ['metrics','peakCall','annotationPlots']});

#--------------------
# PODS
#--------------------
## Here starts the pipeline steps documentation, please change it accordingly any time you add/remove/modify a step

=head1 CHIPSEQ PIPELINE STEPS

The standard peak detection analysis performs the following steps:

B<trimming> :  Raw reads quality trimming and removing of Illumina adapters is performed using trimmomatic.

B<alignment> :  The filtered reads are aligned to a reference genome. The alignment is done per lane of sequencing, and then merged for a complete Binary Alignment Map file (.bam). The alignment software used is bwa , and the merging and sorting is done with the picard software.

B<metrics> :  The number of raw/filtered and aligned reads per sample are computed at this stage.

B<qcTagDirectories> :  The Homer Tag directories, used to check for quality metrics are computed at this step. 

B<qcPlots> :  Sequencing quality metrics as tag count, tag autocorrelation, sequence bias and GC bias are generated. 

B<wiggle> :  Wiggle Track Format files are generated from the aligned reads using Homer. The resulting files can be loaded in browsers like IGV or UCSC.

B<peakCall> :  Peaks are called using the MACS2 software. Different calling strategies are used for narrow and broad peaks. The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run. The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100. The default mfold parameter of MACS is [10,30].  

B<annotation> : The peaks called previously are annotated with HOMER using RefSeq annotations for the reference genome. Gene ontology and genome ontology analysis are also performed at this stage.

B<motif> :  De novo and known motif analysis per design are performed using HOMER .

B<annotationPlots>: The peak location statistics. The following peak location statistics are generated per design: proportions of the genomic locations of the peaks. The locations are: Gene (exon or intron), Proximal ([0;2] kb upstream of a transcription start site), Distal ([2;10] kb upstream of a transcription start site), 5d ([10;100] kb upstream of a transcription start site), Gene desert (>= 100 kb upstream or downstream of a transcription start site), Other (anything not included in the above categories); The distribution of peaks found within exons and introns; The distribution of peak distance relative to the transcription start sites (TSS); the Location of peaks per design.


B<deliverable> : The Standard report. A summary html report is automatically generated by the pipeline. This report contains description of the sequencing experiment as well as a detailed presentation of the pipeline steps and results. Various Quality Control (QC) summary statistics are included in the report and additional QC analysis is accessible for download directly through the report. The report includes also the main references of the software and methods used during the analysis, together with the full list of parameters passed to the pipeline main script.

=cut end of documentation

my %globalDep;
for my $stepName (@steps) {
  $globalDep{$stepName -> {'name'} } ={};
}

my $designFilePath;
my $workDir;
my $currentWorkDir;
my $configFile;


&main();

sub printUsage {
  print "Version: ".$Version::version."\n";
  print "\nUsage: perl ".$0." \n";
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
  getopts('c:s:n:d:w:', \%opts);

  if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'n'}) || !defined($opts{'d'}) || !defined($opts{'w'})) {
    printUsage();
    exit(1);
  }

  my %jobIdVarPrefix;
  my %cfg = LoadConfig->readConfigFile($opts{'c'});
  my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($opts{'n'});
  $designFilePath = $opts{'d'};
  # get Design groups
  my $rHoAoA_designGroup = MACS2::getDesign(\%cfg,$designFilePath);
  $workDir = abs_path($opts{'w'});
  $currentWorkDir = getcwd();
  $configFile =  abs_path($opts{'c'});
  #generate sample jobIdprefix
  my $cpt = 1;

  for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
    my $cpt2 = 1;
    $jobIdVarPrefix{$sampleName} = $cpt;
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
    for my $rH_laneInfo (@$rAoH_sampleLanes) {
      $jobIdVarPrefix{$sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} = $cpt . '_' . $cpt2;
      $cpt2++;
    }
    $cpt++;
  }
  $cpt = 1;
  
  # generate design jobIdprefix
  for my $designName (keys %{$rHoAoA_designGroup}) {
    $jobIdVarPrefix{$designName} = $cpt;
    $cpt++;
  }

  # List user-defined step index range.
  # Shift 1st position to 0 instead of 1
  my @stepRange = map($_ - 1, parse_range($opts{'s'}));

  SubmitToCluster::initPipeline($workDir);

  for my $current (@stepRange) {
    my $fname = $steps[$current]->{'name'};
    my $loopType = $steps[$current]->{'stepLoop'};
    my $subref = \&$fname;
    my @aOjobIDList=();
    
    if ($loopType eq 'sample') {
      for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
        my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
        # Tests for the first step in the list. Used for dependencies.
        my $jobIdVar = &$subref($current , \%cfg, $sampleName, $rAoH_sampleLanes, \%jobIdVarPrefix);
        if (defined($jobIdVar)) {
          $globalDep{$fname}->{$sampleName} = $jobIdVar;
          # This is for global steps depending on sample / design  steps
          if(defined ($globalDep{$fname}->{'experiment'})){
            @aOjobIDList=@{$globalDep{$fname}->{'experiment'}};
          }
          push(@aOjobIDList, $jobIdVar);
          $globalDep{$fname}->{'experiment'} =  \@aOjobIDList;
        }
      }
    } elsif ($loopType eq 'design') {
      for my $design (keys %{$rHoAoA_designGroup}) {        
        my $jobIdVar = &$subref($current, \%cfg, $design, $rHoAoA_designGroup, $rHoAoH_sampleInfo, \%jobIdVarPrefix);
        if (defined($jobIdVar)) {
          $globalDep{$fname}->{$design} = $jobIdVar;
          # This is for global steps depending on sample / design  steps
          if(defined ($globalDep{$fname}->{'experiment'})){
            @aOjobIDList=@{$globalDep{$fname}->{'experiment'}};
          }
          push(@aOjobIDList, $jobIdVar);
          $globalDep{$fname}->{'experiment'} =  \@aOjobIDList;
        }
      }
    } elsif($loopType eq 'experiment') {      
      my $jobIdVar = &$subref($current, \%cfg, $designFilePath, $rHoAoA_designGroup, $rHoAoH_sampleInfo, \%jobIdVarPrefix );
      if (defined($jobIdVar)) {
        if(defined ($globalDep{$fname}->{'experiment'})){
           @aOjobIDList=@{$globalDep{$fname}->{'experiment'}};
        }
        push(@aOjobIDList, $jobIdVar);
        $globalDep{$fname}->{'experiment'} = \@aOjobIDList;
      }
    } else {
      # Tests for the first step in the list. Used for dependencies.
      my $jobIdVar = &$subref($current, \%cfg, $designFilePath, $rHoAoA_designGroup, $rHoAoH_sampleInfo, \%jobIdVarPrefix);
      if (defined($jobIdVar)) {
        $globalDep{$fname}->{$fname} = $jobIdVar;
      }
    }
  }
}

sub bamToFastq {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rH_jobIdPrefixe = shift;  

  # BamToFastq job IDS per sample
  my $bamToFastqJobIdVarNameSample = undef;

  # samples per lane
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
        die "Error in rnaSeqDeNovoAssembly::bamToFastq: unknown run type (can be 'SINGLE_END' or 'PAIRED_END' only): " . $rH_laneInfo->{'runType'};
      }
    }
    if(!$rO_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "bamToFastq", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, "BAMTOFASTQ", $jobDependency, $sampleName, $rO_job);

      if(!defined($bamToFastqJobIdVarNameSample)) {
        $bamToFastqJobIdVarNameSample = $ro_job->getCommandJobId(0);
      } else {
        $bamToFastqJobIdVarNameSample .= LoadConfig::getParam($rH_cfg, 'bamToFastq', 'clusterDependencySep') . $ro_job->getCommandJobId(0);
      }
    }
  }
  return $bamToFastqJobIdVarNameSample;
}

sub trimming {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rH_jobIdPrefixe = shift;  

  # Trimming job IDS per sample
  my $trimJobIdVarNameSample = undef;
  
  # samples per lane
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $outputDir = 'reads/'.$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    print 'mkdir -p '.$outputDir."\n";
    # Run trimmomatic  
    my $ro_job = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo, $outputDir);
    if(!$ro_job->isUp2Date()) {
      SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, "TRIM".$rH_jobIdPrefixe ->{$sampleName.'.' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}} , undef, $sampleName, $ro_job);
      if(!defined($trimJobIdVarNameSample)) {
        $trimJobIdVarNameSample = $ro_job->getCommandJobId(0);
      } else {
        $trimJobIdVarNameSample .= LoadConfig::getParam($rH_cfg, 'trim', 'clusterDependencySep') . $ro_job->getCommandJobId(0);
      }
    }

  }
  return $trimJobIdVarNameSample; 
}


sub aligning{
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rH_jobIdPrefixe = shift;

  my $jobDependency = undef;
  my $setJobId = 0;
  
  # Control dependencies
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if( defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }
  
  print "BWA_JOB_IDS=\"\"\n";
  print "mkdir -p metrics\n";
  
  
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $pair1;
    my $pair2;
    my $single;
    my $rgId = $rH_laneInfo->{'libraryBarcode'} . "_" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rgSampleName = $rH_laneInfo->{'name'};
    my $rgLibrary = $rH_laneInfo->{'libraryBarcode'};
    my $rgPlatformUnit = 'run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    my $rgCenter = LoadConfig::getParam($rH_cfg, 'aln', 'bwaInstitution');
    my $qfilterRead=LoadConfig::getParam($rH_cfg, 'aln', 'filterReadsMAPQ', 0, 'int');
    
    # Threshold for MAPQ to filter reads
    if (!defined($qfilterRead) || $qfilterRead < 1 ){
      $qfilterRead=15;
    }

    # Align lanes
    my $outputAlnDir = 'alignment/' . $sampleName . '/run' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    print 'mkdir -p ' . $outputAlnDir . "\n";
    my $outputAlnPrefix = $outputAlnDir.'/'.$sampleName.'.'.$rH_laneInfo->{'libraryBarcode'};
    my $bwaJobId = undef;

    if ($rH_laneInfo->{'runType'} eq "SINGLE_END") {
      $single = 'reads/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName . '.' .$rH_laneInfo->{'libraryBarcode'}. '.t' . LoadConfig::getParam($rH_cfg, 'trim','minQuality', 1, 'int') . 'l' . LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int') . '.single.fastq.gz';
      $pair1 = undef;
      $pair2 = undef;

    } elsif ($rH_laneInfo->{'runType'} eq "PAIRED_END") {
      $single = undef;
      $pair1 = 'reads/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName . '.' .$rH_laneInfo->{'libraryBarcode'}. '.t' . LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int') . 'l' . LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int') . '.pair1.fastq.gz';
      $pair2 = 'reads/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName . '.' .$rH_laneInfo->{'libraryBarcode'}. '.t' . LoadConfig::getParam($rH_cfg, 'trim', 'minQuality', 1, 'int') . 'l' . LoadConfig::getParam($rH_cfg, 'trim', 'minLength', 1, 'int') . '.pair2.fastq.gz';
    }
    # BWA 
    my $ro_bwaJob = BWA::aln($rH_cfg, $sampleName, $pair1, $pair2, $single, $outputAlnPrefix, $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter);
    if(!$ro_bwaJob->isUp2Date()) {
      if($ro_bwaJob->getNbCommands() == 3) {        
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read1.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ1ALN', $jobDependency, $sampleName, $ro_bwaJob, 0 );
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read2.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ2ALN', $jobDependency, $sampleName, $ro_bwaJob, 1 );
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'sampe.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA', $ro_bwaJob->getCommandJobId(0).LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$ro_bwaJob->getCommandJobId(1), $sampleName, $ro_bwaJob, 2 );
        print 'BWA_JOB_IDS=${BWA_JOB_IDS#:}'.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$ro_bwaJob->getCommandJobId(2)."\n";
        $bwaJobId = $ro_bwaJob->getCommandJobId(2);
        $setJobId = 1;
      } else {
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READALN', $jobDependency, $sampleName, $ro_bwaJob, 0 );
        SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'samse.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA',  $ro_bwaJob->getCommandJobId(0), $sampleName, $ro_bwaJob, 1 );
        print 'BWA_JOB_IDS=${BWA_JOB_IDS#:}'.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$ro_bwaJob->getCommandJobId(1)."\n";
        $bwaJobId = $ro_bwaJob->getCommandJobId(1);
        $setJobId = 1;
      }
    }
    # Filter uniquely mapped reads
    my $dir = 'alignment/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    my $InputBamFile = $dir . $rH_laneInfo->{'name'}. '.' .$rH_laneInfo->{'libraryBarcode'}. ".sorted.bam";
    my $OutputBamFile = $dir . $rH_laneInfo->{'name'} . '.' .$rH_laneInfo->{'libraryBarcode'}. ".sorted.filtered.bam";
    my $ro_job = SAMtools::viewFilter($rH_cfg, $InputBamFile, " -b -F4 -q " . $qfilterRead . " ", $OutputBamFile);
    if(!$ro_job->isUp2Date()) {
     SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'filter.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'SAMTOOLSFILTER',  $bwaJobId , $sampleName, $ro_job);
     print 'BWA_JOB_IDS=${BWA_JOB_IDS#:}' . LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep') . $ro_job->getCommandJobId(0) . "\n\n";
    }

  }
 
  # Merge /sort reads
  if($setJobId ==0){
    $jobDependency = undef;
  }else{
    $jobDependency = '$BWA_JOB_IDS';
  }
  my $mergeJobId = undef ;
  my @inputBams;
  my $outputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.bam';
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $dir = 'alignment/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    my $sortedLaneBamFile =$dir . $rH_laneInfo->{'name'} . '.' .$rH_laneInfo->{'libraryBarcode'}. ".sorted.filtered.bam" ;
    push(@inputBams, $sortedLaneBamFile);
  }
  my $ro_job  = Picard::mergeFiles($rH_cfg, $sampleName, \@inputBams, $outputBAM);
  if(!$ro_job->isUp2Date()) {  
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeLanes", undef, 'MERGELANES' . $rH_jobIdPrefixe->{$sampleName}, $jobDependency, $sampleName, $ro_job);
    $mergeJobId = $ro_job->getCommandJobId(0);
  }
 
  # Mark as duplicates
  my $markDupJobID=undef;
  my $mergedBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.merged.bam';
  $outputBAM = 'alignment/' . $sampleName . '/' . $sampleName . '.sorted.bam';
  my $outputMetrics = 'alignment/' . $sampleName . '/' . $sampleName . '.sorted.metrics';
  my $rO_job = Picard::markDup($rH_cfg, $sampleName, $mergedBAM , $outputBAM , $outputMetrics);
  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'SAMPLEMARKDUP'. $rH_jobIdPrefixe->{$sampleName}, $mergeJobId , $sampleName, $rO_job);
    $markDupJobID=$rO_job->getCommandJobId(0);
  }
  return $markDupJobID;
}

sub metrics {
  my $stepId              = shift;
  my $rH_cfg              = shift;
  my $sampleName          = shift;
  my $rAoH_sampleLanes    = shift;
  my $rH_jobIdPrefixe     = shift;

  my $jobDependency       = undef;
  my $flagstatJobId       = undef;
  my $metricsJobIDs  = undef;

  # Control dependencies
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }
  
  # Generate alignment stats
  for my $rH_laneInfo (@$rAoH_sampleLanes) {

   # Get trimmed read count is automatically done by trimmomatic
    my $groupName       = $sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}.'.' .$rH_laneInfo->{'libraryBarcode'};
    my $trimStatsFile   =  'reads/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/" . $sampleName . '.' .$rH_laneInfo->{'libraryBarcode'}. '.trim.stats.csv';

    # Compute flagstats per sample / lane
    my $dir = 'alignment/' . $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    my $inputBamFile = $dir . $sampleName . '.' .$rH_laneInfo->{'libraryBarcode'}. '.sorted.bam';
    my $flagStatsFile = $dir . $sampleName . '.' .$rH_laneInfo->{'libraryBarcode'} . '.sorted.bam.flagstat';

    my $ro_job  = SAMtools::flagstat($rH_cfg, $inputBamFile, $flagStatsFile);
    if(!$ro_job->isUp2Date()) {  
      SubmitToCluster::printSubmitCmd($rH_cfg, "flagstat", undef, 'FLAGSTAT', $jobDependency, $groupName, $ro_job);
      $flagstatJobId = $ro_job->getCommandJobId(0);
    }

    ## Merge read stats
    my $outputFile= 'metrics/' . $groupName . '.readstats.csv';
    my $sampleColumnName       = $sampleName . ' ' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}.' ' .$rH_laneInfo->{'libraryBarcode'};
    $ro_job  = Metrics::mergePrintReadStats($rH_cfg, $sampleColumnName, $trimStatsFile, $flagStatsFile, $outputFile);
    if(!$ro_job->isUp2Date()) {  
      SubmitToCluster::printSubmitCmd($rH_cfg, "mergeMetrics", undef, 'MERGEREADSTATLANE' . $rH_jobIdPrefixe->{$sampleName . '.' . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}}, $flagstatJobId, $groupName, $ro_job);
      if(defined($metricsJobIDs)){
        $metricsJobIDs=$metricsJobIDs.LoadConfig::getParam($rH_cfg, 'metrics', 'clusterDependencySep').$ro_job->getCommandJobId(0);    
      }else{
        $metricsJobIDs=$ro_job->getCommandJobId(0);    
      }  
    }
  }
  # Run general metrics for merged /marked as duplicate bam files
  my $dir = 'alignment/' . $sampleName . "/";
  my $inputBamFile = $dir . $sampleName . '.sorted.bam';
  my $flagStatsFile = $dir . $sampleName . '.sorted.bam.flagstat';

  my $ro_job  = SAMtools::flagstat($rH_cfg, $inputBamFile, $flagStatsFile);
  if(!$ro_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "flagstat", undef, 'FLAGSTAT', $jobDependency, $sampleName, $ro_job);
    $flagstatJobId = $ro_job->getCommandJobId(0);
  }
  ## Parse flagstats
  my $outputFile= 'metrics/' . $sampleName . '.alnstats.csv';    
  $ro_job  = Metrics::mergePrintReadStats($rH_cfg, $sampleName , "", $flagStatsFile , $outputFile);
  if(!$ro_job->isUp2Date()) {  
    SubmitToCluster::printSubmitCmd($rH_cfg, "mergeMetrics", undef, 'MERGEREADSTATSAMPLE' . $rH_jobIdPrefixe->{$sampleName}, $flagstatJobId, $sampleName,  $ro_job);
    if(defined($metricsJobIDs)){
        $metricsJobIDs=$metricsJobIDs.LoadConfig::getParam($rH_cfg, 'metrics', 'clusterDependencySep').$ro_job->getCommandJobId(0);    
    }else{
        $metricsJobIDs=$ro_job->getCommandJobId(0);    
    }  
  }
  return $metricsJobIDs;
}

sub qcTagDirectories {
  my $stepId      = shift;
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rAoH_sampleLanes = shift;
  my $rH_jobIdPrefixe = shift;

  my $jobDependency = undef;
  my $qcTagsJobID   = undef;

  # Control dependencies
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }
  
  
  my $inputBAM = 'alignment/' . $sampleName . "/" . $sampleName . ".sorted.bam";
  # Create command
  
  my $ro_job = HOMER::makeTagDirectory($rH_cfg, $sampleName, $inputBAM, 'tags');
  if(!$ro_job->isUp2Date()) { 
    SubmitToCluster::printSubmitCmd($rH_cfg, "qcTags", undef, 'QCTAGS' . $rH_jobIdPrefixe->{$sampleName}, $jobDependency, $sampleName, $ro_job );
    $qcTagsJobID = $ro_job->getCommandJobId(0);
  }
  
  return $qcTagsJobID;
}

sub wiggle {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes = shift;
  my $rH_jobIdPrefixe = shift;

  my $jobDependency = undef;

  # Control dependencies
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$sampleName})) {
    $jobDependency = $globalDep{$parentStep}->{$sampleName};
  }

  my $tagDir         =  'tags/' . $sampleName;
  my $outputWiggle   =  'tracks/' . $sampleName . '/' . $sampleName . '.ucsc.bedGraph.gz';
  my $prefixJobName  =  undef;

  print "mkdir -p tracks/$sampleName\n";

  my $wiggleJobId = undef;
  my $ro_job  = HOMER::makeUCSCFile($rH_cfg, $sampleName, $tagDir, $outputWiggle);
  if(!$ro_job->isUp2Date()) { 
    SubmitToCluster::printSubmitCmd($rH_cfg, "wiggle", $prefixJobName, 'WIGGLE' . $rH_jobIdPrefixe->{$sampleName}, $jobDependency, $sampleName, $ro_job);
    $wiggleJobId  = $ro_job->getCommandJobId(0);
  }
  return $wiggleJobId;
}

sub qcPlots {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $designFilePath = shift;

  my $jobDependency = undef;
  
  # Control dependencies
  my $parentStep = $steps[$stepId]->{'parentStep'};
  my $jobDependencies = "";
  foreach my $stepName (@{$parentStep}) {
    if(defined($globalDep{$stepName}->{'experiment'})){
      my @aOjobIDList=@{$globalDep{$stepName}->{'experiment'}};
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), @aOjobIDList);
    }
  }
  if(length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }
  $jobDependency=$jobDependencies;

  print "mkdir -p graphs\n";
  my $qcPlotsJobId = undef;

  my $ro_job  = HOMER::qcPlotsR($rH_cfg, $designFilePath, $workDir);
  if(!$ro_job->isUp2Date()) { 
    SubmitToCluster::printSubmitCmd($rH_cfg, "qcGraphs", undef, 'QCPLOTS', $jobDependency, undef, $ro_job);
    $qcPlotsJobId = $ro_job->getCommandJobId(0);
  }
  return $qcPlotsJobId;
}

sub peakCall {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $design = shift;
  my $rHoAoA_designGroup = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rH_jobIdPrefixe = shift;

  my $jobDependency = undef;
  my $PeakCallJobIdVarNameSample = undef;
  my $parentStep = $steps[$stepId]->{'parentStep'};

  ## Design
  # Process Treatments (and controls when applicable)
  my $numberTreatments =  @{$rHoAoA_designGroup->{$design}->[0]};
  my $numberControls   =  @{$rHoAoA_designGroup->{$design}->[1]};
  my $control;
  my $treatment;
  my $controlBam;
  my $treatmentBam;

  if ($numberTreatments >= 1) {
    # At least one treatment
    for (my $j = 0; $j < $numberTreatments; $j++) {
      if ($numberControls == 1) {
          $control       = $rHoAoA_designGroup->{$design}->[1]->[0];
          $controlBam    = 'alignment/' . $control . '/' . $control . '.sorted.bam';
      } elsif ($numberControls == $numberTreatments) {
          $control       = $rHoAoA_designGroup->{$design}->[1]->[$j];
          $controlBam    = 'alignment/' . $control . '/' .$control . '.sorted.bam';
      } else {
          $control       = undef;
          $controlBam    = undef;
      }
      $treatment    =  $rHoAoA_designGroup->{$design}->[0]->[$j];
      $treatmentBam =  'alignment/' . $treatment . '/' . $treatment . '.sorted.bam';
      my $outputPath ;        
      if (defined($control)) {
        print "# design " . $design . ", treatment= " . $treatment . ", control=" . $control . "\n";
        $outputPath= 'peak_call/' . $design. '__' . $control . '__' . $treatment;
      } else {
        print "# design " . $design . ", treatment= " . $treatment . "" . "\n";
        $outputPath= 'peak_call/' . $design. '__' . '' . '__' . $treatment;
      }
      # Run type for treatment
      my $paired = undef;
      my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$treatment};
      for my $rH_laneInfo (@$rAoH_sampleLanes) {
        if ($rH_laneInfo->{'runType'} eq "SINGLE_END") {
          $paired = 0;
        } elsif ($rH_laneInfo->{'runType'} eq "PAIRED_END") {
          $paired = 1;
        } else {
          die "Unknown runType: " . $rH_laneInfo->{'runType'} . "\n";
        }
        last if (defined($paired));
      }
      if (!defined($rHoAoH_sampleInfo->{$treatment})) {
        die "ERROR: Sample " . $treatment . " in design file is not present in NANUQ sample file" . "\n";
      } elsif (defined($control) && !defined($rHoAoH_sampleInfo->{$control})) {
        die "ERROR: Sample " . $control . " in design file is not present in NANUQ sample file " . "\n";
      }
      # MACS command
      my $ro_job  = MACS2::generatePeaks($rH_cfg, $design, $treatmentBam, $controlBam, $rHoAoA_designGroup->{$design}->[2]->[0], $outputPath, $paired);
      if(!$ro_job->isUp2Date()) { 
        if (defined($control) && defined($globalDep{$parentStep}->{$control}) && defined($globalDep{$parentStep}->{$treatment})) {
          $jobDependency = $globalDep{$parentStep}->{$control} . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep') . $globalDep{$parentStep}->{$treatment};
        } elsif ( defined($globalDep{$parentStep}->{$treatment})) {
          $jobDependency = $globalDep{$parentStep}->{$treatment};
        }
        print "mkdir -p $outputPath\n";
        SubmitToCluster::printSubmitCmd($rH_cfg, 'peakCall', $j, 'MACS2PEAKCALL' . $rH_jobIdPrefixe->{$design} . $j, $jobDependency, $design, $ro_job);
        $PeakCallJobIdVarNameSample .= $ro_job->getCommandJobId(0) . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
      }
    }
  } else {
    ; # do nothing
  }
  if (defined($PeakCallJobIdVarNameSample)) {
    $PeakCallJobIdVarNameSample = substr $PeakCallJobIdVarNameSample, 0, -1;
  }
  return $PeakCallJobIdVarNameSample;
}

sub annotation {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $design = shift;
  my $rHoAoA_designGroup = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rH_jobIdPrefixe = shift;

  my $jobDependency = undef;
  my $setJobId = 0;
  # Control dependencies
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$design})) {
    $jobDependency = $globalDep{$parentStep}->{$design};
  }
  ##Iterate over design
  my $numberTreatments =  @{$rHoAoA_designGroup->{$design}->[0]};
  my $numberControls   =  @{$rHoAoA_designGroup->{$design}->[1]};
  my $control;
  my $treatment;

  my $annotationJobIdGroups = undef;
  if ($numberTreatments >= 1) {
    print "mkdir -p annotation\n";
    # At least one treatment
    for (my $j = 0; $j < $numberTreatments; $j++) {
      if ($numberControls == 1) {
          $control       = $rHoAoA_designGroup->{$design}->[1]->[0];
      } elsif ($numberControls == $numberTreatments) {
          $control       = $rHoAoA_designGroup->{$design}->[1]->[$j];
      } else {
          $control       = "";
      }
      $treatment    =  $rHoAoA_designGroup->{$design}->[0]->[$j];
      my $peakCallPath = 'peak_call/' . $design. '__' . $control . '__' . $treatment;
      my $outputPath = 'annotation/' . $design. '__' . $control . '__' . $treatment;
      my $peakType   =  $rHoAoA_designGroup->{$design}->[2]->[0];
      my $annotationJobID = undef;
      if ($peakType eq "N") {
        my $bedName = $peakCallPath . '/' . $design . '_peaks.bed';
        my $ro_job  =  HOMER::annotatePeaks($rH_cfg, $design, $bedName, $outputPath);
        if(!$ro_job->isUp2Date()) { 
          print "mkdir -p $outputPath\n";
          SubmitToCluster::printSubmitCmd($rH_cfg, 'annotation', undef, 'HOMERANNOTATION' . $rH_jobIdPrefixe->{$design}, $jobDependency, $design, $ro_job);          
          $annotationJobID=$ro_job->getCommandJobId(0);
          $setJobId = 1;
        }
        # Create annotation stats
        $ro_job  = HOMER::annotateStats($rH_cfg, $design, $outputPath. '/' . $design . '.annotated.csv' , $outputPath. '/' . $design);
        if(!$ro_job->isUp2Date()) { 
            SubmitToCluster::printSubmitCmd($rH_cfg, 'annotationStats', undef, 'HOMERANNOTATIONSTATS' . $rH_jobIdPrefixe->{$design}, $annotationJobID , $design, $ro_job);          
            if (!defined($annotationJobIdGroups)){
              $annotationJobIdGroups = $ro_job->getCommandJobId(0);
            }else{
              $annotationJobIdGroups .= LoadConfig::getParam($rH_cfg, 'annotation', 'clusterDependencySep'). $ro_job->getCommandJobId(0);
            }
            $setJobId = 1;
          }
      } else {
        ; # do nothing
      }
    }
  }
  if($setJobId ==0){
    return undef;
  }
  return $annotationJobIdGroups;
}

sub motif {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $design = shift;
  my $rHoAoA_designGroup = shift;
  my $rHoAoH_sampleInfo = shift;
  my $rH_jobIdPrefixe = shift;

  my $jobDependency = undef;

  # Control dependencies
  my $parentStep = $steps[$stepId]->{'parentStep'};
  if(defined($globalDep{$parentStep}->{$design})) {
    $jobDependency = $globalDep{$parentStep}->{$design};
  }
  # Iterate over design
  my $numberTreatments =  @{$rHoAoA_designGroup->{$design}->[0]};
  my $numberControls   =  @{$rHoAoA_designGroup->{$design}->[1]};
  my $control;
  my $treatment;
  
  my $motifJobIdGroups = undef;
  if ($numberTreatments >= 1) {
    # At least one treatment
    for (my $j = 0; $j < $numberTreatments; $j++) {
      if ($numberControls == 1) {
          $control       = $rHoAoA_designGroup->{$design}->[1]->[0];
      } elsif ($numberControls == $numberTreatments) {
          $control       = $rHoAoA_designGroup->{$design}->[1]->[$j];
      } else {
          $control       = "";
      }
      $treatment    =  $rHoAoA_designGroup->{$design}->[0]->[$j];
      my $peakCallPath = 'peak_call/' . $design. '__' . $control . '__' . $treatment;
      my $outputPath = 'annotation/' . $design. '__' . $control . '__' . $treatment. '/' . $design;
      my $peakType =  $rHoAoA_designGroup->{$design}->[2]->[0];
      if ($peakType eq "N") {
        my $bedName = $peakCallPath . '/' . $design . '_peaks.bed';
        my $ro_job  =  HOMER::generateMotif($rH_cfg, $design, $bedName, $outputPath);
        if(!$ro_job->isUp2Date()) { 
          print "mkdir -p $outputPath\n";
          SubmitToCluster::printSubmitCmd($rH_cfg, 'motif', undef, 'HOMERMOTIF' . $rH_jobIdPrefixe->{$design}, $jobDependency, $design, $ro_job);          
          $motifJobIdGroups .= $ro_job->getCommandJobId(0) . LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep');
        }
      } else {
        ; # do nothing
      }
    }
  }
  return $motifJobIdGroups;
}

sub annotationPlots {
  my $stepId = shift;
  my $rH_cfg = shift;
  my $designFilePath = shift;
  my $jobDependency = undef;
  
  # Control dependencies
  my $parentStep = $steps[$stepId]->{'parentStep'};
  my $jobDependencies = "";
  foreach my $stepName (@{$parentStep}) {
    if(defined($globalDep{$stepName}->{'experiment'})){
      my @aOjobIDList=@{$globalDep{$stepName}->{'experiment'}};
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'). join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), @aOjobIDList);
    }
  }
  if(length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }  
  $jobDependency=$jobDependencies;
  
  print "mkdir -p graphs\n";
  my $annotPlotsJobId = undef;

  my $ro_job  = HOMER::annotatePlotsR($rH_cfg, $designFilePath, $workDir);
  if(!$ro_job->isUp2Date()) { 
    SubmitToCluster::printSubmitCmd($rH_cfg, "annotationPlots", undef, 'ANNOTATIONPLOTS', $jobDependency, undef, $ro_job);
    $annotPlotsJobId = $ro_job->getCommandJobId(0);
  }
  return $annotPlotsJobId;

  return 1;
}

sub deliverable{
  my $stepId = shift;
  my $rH_cfg = shift;
  
  my $parentStep = $steps[$stepId]->{'parentStep'};
  my $jobDependencies = "";
  foreach my $stepName (@{$parentStep}) {
    if(defined($globalDep{$stepName}->{'experiment'})){
      my @aOjobIDList=@{$globalDep{$stepName}->{'experiment'}};
      $jobDependencies .= LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'). join(LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), @aOjobIDList);
    }
  }
  if(length($jobDependencies) == 0) {
    $jobDependencies = undef;
  } else {
    $jobDependencies = substr($jobDependencies, 1);
  }

  my $rO_job = GqSeqUtils::clientReport($rH_cfg,  $configFile, $currentWorkDir, 'CHIPseq') ;

  if(!$rO_job->isUp2Date()) {
    SubmitToCluster::printSubmitCmd($rH_cfg, "deliverable", 'REPORT', 'CHIPSEQREPORT', $jobDependencies , 'allSamples', $rO_job);
  }

  return $rO_job->getCommandJobId(0);
}

1;
