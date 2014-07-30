#!/usr/bin/env perl

=head1 NAME

I<illuminaRunProcessing>

=head1 SYNOPSIS

illuminaRunProcessing.pl

=head1 DESCRIPTION

B<illuminaRunProcessing> Pipeline to process raw data from Illumina NGS sequencers.

=head1 AUTHOR

B<Marc Michaud> - I<marc.michaud@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Parse::Range> Step range parsing

B<XML::Simple> for RunInfo.xml file parsing

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
use XML::Simple;
use Data::Dumper;
use Parse::Range qw(parse_range);
use File::Basename;

use CountIlluminaBarcodes;
use BVATools;
use BWA;
use LoadConfig;
use Picard;
use SampleSheet;
use SubmitToCluster;
use Tools;
use Version;
use Metrics;


my @steps;
push(@steps, {'name' => 'generateIndexCount', 'stepLoop' => 'lane', 'parentStep' => []});
push(@steps, {'name' => 'generateFastq', 'stepLoop' => 'lane', 'parentStep' => []});
push(@steps, {'name' => 'generateMD5', 'stepLoop' => 'sample', 'parentStep' => ['generateFastq']});
push(@steps, {'name' => 'generateQCGraphs', 'stepLoop' => 'sample', 'parentStep' => ['generateFastq']});
push(@steps, {'name' => 'generateBlasts', 'stepLoop' => 'sample', 'parentStep' => ['generateFastq']});
push(@steps, {'name' => 'align', 'stepLoop' => 'sample', 'parentStep' => ['generateFastq']});
push(@steps, {'name' => 'laneMetrics', 'stepLoop' => 'sample', 'parentStep' => ['align']});
push(@steps, {'name' => 'generateBamMd5', 'stepLoop' => 'sample', 'parentStep' => ['align']});
push(@steps, {'name' => 'startCopyNotification' , 'stepLoop' => 'lane' , 'parentStep' => ['generateIndexCount','generateMD5','generateQCGraphs','generateBlasts','laneMetrics','generateBamMd5']});
push(@steps, {'name' => 'copy' , 'stepLoop' => 'lane' , 'parentStep' => ['generateIndexCount','generateMD5','generateQCGraphs','generateBlasts','laneMetrics','generateBamMd5']});
push(@steps, {'name' => 'endCopyNotification' , 'stepLoop' => 'lane' , 'parentStep' => ['copy']});

my $UNALIGNED_DIR="Unaligned";
my $ALIGNED_DIR="Aligned";

# Create step hash indexed by step name for easy retrieval
my %H_steps =  map {$_->{'name'} => $_} @steps;

# Global scope variables
my $GLOBAL_DEP_KEY = "#global";
my %generatedIntervalLists;

my $rHoH_genomes;
my $rHoH_defaultGenomes;
my $casavaSheet;
my $globalNumberOfMismatch;
my $mask;
my $firstIndex;
my $lastIndex;

&main();

sub printUsage {
  print "Version: ".$Version::version."\n";
  print "\nUsage: perl ".$0." -c config.ini -s steps -l nb -r /path/to/run\n";
  print "\t-c  config file\n";
  print "\t-s  step range e.g. '1,3', '2-5', '1,4-7,11'\n";  
  print "\t-l  lane number\n";
  print "\t-r  run directory\n";
  print "\t-n  nanuq sample sheet. Optional, default=RUNDIRECTORY/run.nanuq.csv\n";
  print "\t-i  Illumina (Casava) sample sheet. Optional, default=RUNDIRECTORY/SampleSheet.nanuq.csv\n";
  print "\t-f  Force download of samples sheets, even if already existing, Optional, default= no force\n";
  print "\t-m  Number of mismatches. Optional, default=1\n";
  print "\t-x  First index nucleotide to use. Optional, default=1\n";
  print "\t-y  Last index nucleotide to use. Optional, default=last\n";
  print "\n";
  print "Steps:\n";
  for(my $idx=0; $idx < @steps; $idx++) {
    print "".($idx+1).'- '.$steps[$idx]->{'name'}."\n";
  }
  print "\n";
  return;
}

sub main {
  my %opts;
  getopts('fc:s:m:n:l:r:i:x:y:', \%opts);

  if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'l'}) || !defined($opts{'r'})) {
    printUsage();
    exit(1);
  }

  my $runDirectory = $opts{'r'};
  my $lane = $opts{'l'};
  my %cfg = LoadConfig->readConfigFile($opts{'c'});
  $globalNumberOfMismatch = defined($opts{'m'}) ? $opts{'m'} : LoadConfig::getParam(\%cfg, 'default', 'numberMismatches');
  $casavaSheet = defined($opts{'i'}) ? $opts{'i'} : $runDirectory . "/SampleSheet.nanuq.csv";
  my $nanuqSheet = defined($opts{'n'}) ? $opts{'n'} : $runDirectory . "/run.nanuq.csv";
  $firstIndex = $opts{'x'};
  $lastIndex = $opts{'y'};
  my $forceDownload = $opts{'f'};
  
  if (!defined($lastIndex) ||  !($lastIndex > 0)) {
    $lastIndex = 999;
  }
  
  if (!defined($firstIndex) || !($firstIndex > 0)) {
    $firstIndex = 1;
  }
  
  my $runID;
  my $runName;
  if ( $runDirectory =~ /.*_\d+HS\d\d[AB]/ ) {
    ($runName,$runID) = $runDirectory =~ /.*\/(\d+_[^_]+_\d+_[^_]+_(\d+)HS.+)/;
  }
  elsif ( $runDirectory =~ /.*\d+_[^_]+_\d+_.+/ ) {
    ($runName,$runID) = $runDirectory =~ /.*\/(\d+_([^_]+_\d+)_.*)/;
  }

  if (!defined($opts{'i'}) || defined($forceDownload)) {
    # Download casava sheet
    if ( (!-e $casavaSheet) || defined($forceDownload)) {
      my $command = formatCommand("config" => \%cfg, "command" => LoadConfig::getParam(\%cfg, 'default', 'fetchCasavaSheetCommand'), "runDirectory" => $runDirectory, "runID" => $runID, "lane" => $lane, "filename" => $casavaSheet);
      system($command);
    }
  }
  
  if (!defined($opts{'n'}) || defined($forceDownload)) {
    # Download nanuq run sheet
    if ( (!-e $nanuqSheet) || defined($forceDownload)) {
      my $command = formatCommand("config" => \%cfg, "command" => LoadConfig::getParam(\%cfg, 'default', 'fetchNanuqSheetCommand'), "runDirectory" => $runDirectory, "runID" => $runID, "lane" => $lane, "filename" => $nanuqSheet);
      system($command);
    }
  }

  
  $UNALIGNED_DIR = LoadConfig::getParam(\%cfg, 'default', 'unalignedDirPrefix');
  $ALIGNED_DIR  = LoadConfig::getParam(\%cfg, 'default', 'alignedDirPrefix');
  
  my ($nbReads, $rAoH_readsInfo) = parseRunInfoFile($runDirectory ."/RunInfo.xml" );

  $mask = getMask($lane, $rAoH_readsInfo, $firstIndex, $lastIndex);
  my $rAoH_samples = generateIlluminaLaneSampleSheet($lane, $runDirectory, $rAoH_readsInfo, $mask);
  my $rHoAoH_infos = SampleSheet::parseSampleSheetAsHashByProcessingId($nanuqSheet,1);

  $rHoH_genomes = getGenomeList(LoadConfig::getParam(\%cfg, 'default', 'genomesHome'));
  $rHoH_defaultGenomes = getDefaultGenomes(LoadConfig::getParam(\%cfg, 'default', 'defaultSpeciesGenome'));

  # Merge informations from the nanuq run sheet with those from the casava sample sheet
  for my $rH_sample (@$rAoH_samples) {
    for my $rH_sampleFromNanuqSheet (@{$rHoAoH_infos->{$rH_sample->{'processingSheetId'}}}) {
      if ($rH_sampleFromNanuqSheet->{'lane'} == $lane) {
        for my $key (keys %$rH_sampleFromNanuqSheet) {
          $rH_sample->{$key} = $rH_sampleFromNanuqSheet->{$key};
        }
      }
    }
  }

  
  
  # List user-defined step index range.
  # Shift 1st position to 0 instead of 1
  my @stepRange = map($_ - 1, parse_range($opts{'s'}));

  SubmitToCluster::initPipeline($runDirectory);

  my $startStep = $opts{'s'};
  my $endStep = $stepRange[$#stepRange];
  # Go through steps and create global or sample jobs accordingly
  for my $currentStep (@stepRange) {
    my $step = $steps[$currentStep];
    my $stepName = $step->{'name'};
    my $stepRef = \&$stepName;
    $step->{'jobIds'}->{$GLOBAL_DEP_KEY} = ();

    &$stepRef($step, \%cfg, $runDirectory, $runID, $lane, $rAoH_readsInfo, $nbReads, $rAoH_samples);
  }

  my $jobIds = join (LoadConfig::getParam(\%cfg, 'default', 'clusterDependencySep'), map {"\$" . $_} @{$steps[$endStep]->{'jobIds'}->{$GLOBAL_DEP_KEY}});
  print 'export FINAL_STEP_JOB_IDS='.$jobIds."\n";
  return;
}

sub generateIndexCount {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;

  my $indexLength = computeIndexLength($rAoH_readsInfo);
  my $baseCallDir = $runDirectory. "/" . LoadConfig::getParam($rH_cfg, 'generateFastq', 'baseCallDir');
  my $runName = basename($runDirectory);

  my $mask = "";
  my $indexPrinted=0;
  for my $rH_readInfo (@{$rAoH_readsInfo}) {
    if($rH_readInfo->{'isIndexed'} eq 'Y') {
      if($indexPrinted == 0) {
        $mask .= $indexLength.'B';
        $indexPrinted=1;
      }
    } elsif ($indexPrinted == 1) {
        last; # Don't write the last read, it saves some time!
    } else {
      $mask .= $rH_readInfo->{'nbCycles'}.'T';
    }
  }

  if( $indexPrinted == 0) {
    print  "# No Indexes, *NOT* Generating index counts\n";
  } else {
    my $jobDependency = undef;
    my $ro_job = CountIlluminaBarcodes::count($rH_cfg, $baseCallDir, $lane, $mask, $runDirectory ."/RunInfo.xml", $runDirectory.'/'.$runName.'_'.$lane.'.metrics', $globalNumberOfMismatch);
    
    if (!$ro_job->isUp2Date()) {
      print  "# Generating index counts\n";    
      my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "generateIndexCount", "$runID.$lane", 'idx_', $jobDependency, undef, $ro_job);
      push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
    }
  }
  return;
}


sub generateFastq {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;

  my $jobDependency = undef;
  my @outputs;
  
  for my $rH_sample (@$rAoH_sample) {
    push(@outputs, getFastqFilename($runDirectory, $lane, $rH_sample, 1));
    if ($nbReads > 1) {
      push(@outputs, getFastqFilename($runDirectory, $lane, $rH_sample, 2));
    }
  }
  
  my $ro_job = new Job();
  $ro_job->testInputOutputs([$casavaSheet],\@outputs);
  if (!$ro_job->isUp2Date()) {
  
    my $casavaCmd = LoadConfig::getParam($rH_cfg, 'generateFastq', 'casavaCmd');
    my $baseCallDir = $runDirectory. "/" . LoadConfig::getParam($rH_cfg, 'generateFastq', 'baseCallDir');
    my $casavaSampleSheetPrefix = LoadConfig::getParam($rH_cfg, 'generateFastq', 'casavaSampleSheetPrefix');

    validateBarcodes($globalNumberOfMismatch, $rAoH_sample);


    print "# Generate Unaligned directory\n";
    print "module load " . LoadConfig::getParam($rH_cfg, 'generateFastq', 'moduleVersion.casava') . "\n";
    if($mask !~ /I/) {
      print  "$casavaCmd --input-dir $baseCallDir --output-dir $runDirectory/$UNALIGNED_DIR.$lane --tiles s_$lane --sample-sheet $runDirectory/$casavaSampleSheetPrefix$lane.csv --fastq-cluster-count 0\n";
    } else {
      print  "$casavaCmd --input-dir $baseCallDir --output-dir $runDirectory/$UNALIGNED_DIR.$lane --tiles s_$lane --sample-sheet $runDirectory/$casavaSampleSheetPrefix$lane.csv --fastq-cluster-count 0 --mismatches $globalNumberOfMismatch --use-bases-mask $mask\n";
    }
    if (LoadConfig::getParam($rH_cfg, 'generateFastq', 'sendNotification')) {
      print formatCommand("config" => $rH_cfg, "command" => LoadConfig::getParam($rH_cfg, 'generateFastq', 'startNotificationCommand'), "runDirectory" => $runDirectory, "runID" => $runID, "lane" => $lane, "mask" => $mask, "mismatches" => $globalNumberOfMismatch) . "\n";
      print "\n";
    }

    
    my $command = 'cd ' . $UNALIGNED_DIR. '.' . $lane . ' && make -j '.LoadConfig::getParam($rH_cfg, 'generateFastq', 'nbThreads');
    $ro_job->addCommand($command);

    my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "generateFastq", "$runID.$lane", 'fastq.', $jobDependency, undef, $ro_job);
    push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
  }
  return;
}

sub generateMD5 {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;

  my $dependencies = getDependencies($step, $rH_cfg);

  for my $rH_sample (@$rAoH_sample) {
    my $ro_job = new Job();
    
    my $read1File = getFastqFilename($runDirectory, $lane, $rH_sample, 1);
    my $read2File = getFastqFilename($runDirectory, $lane, $rH_sample, 2);
    
    if($nbReads == 1) {
      $ro_job->testInputOutputs([$read1File],[$read1File.'.md5']);
    } else {
      $ro_job->testInputOutputs([$read1File, $read2File],[$read1File.'.md5', $read2File.'.md5']);
    }

    if (!$ro_job->isUp2Date()) {
      my $command = 'md5sum -b '. $read1File . ' > '. $read1File. '.md5';

      if($nbReads > 1) {
        $command .= '; md5sum -b '. $read2File . ' > '. $read2File . '.md5';
      }

      $ro_job->addCommand($command);
      my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "generateMD5", "$runID.$lane", 'md5_'.$rH_sample->{'processingSheetId'}, $dependencies, $rH_sample->{'processingSheetId'}, $ro_job);
      push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
    }
  }
  return;
}

sub generateQCGraphs {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;


  my $dependencies = getDependencies($step, $rH_cfg);

  for my $rH_sample (@$rAoH_sample) {
    my $output = $runDirectory.'/' . $UNALIGNED_DIR . '.'.$lane.'/Project_nanuq/Sample_'.$rH_sample->{'processingSheetId'}.'/qc';
    my $regionName = $rH_sample->{'processingSheetId'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'};
    my $read1 = getFastqFilename($runDirectory, $lane, $rH_sample, 1);
    my $read2;
    if($nbReads > 1) {
      $read2 = getFastqFilename($runDirectory, $lane, $rH_sample, 2);
    }
    my $ro_job = BVATools::qc($rH_cfg, $read1, $read2, "FASTQ", $regionName, $output);
    if (!$ro_job->isUp2Date()) {
      print 'mkdir -p '.$runDirectory.'/' . $UNALIGNED_DIR. '.'.$lane.'/Project_nanuq/Sample_'.$rH_sample->{'processingSheetId'}."/qc\n";
      my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "generateQCGraphs", "$runID.$lane", 'qc_'.$rH_sample->{'processingSheetId'}, $dependencies, $rH_sample->{'processingSheetId'}, $ro_job);
      push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
    }
  }
  return;
}

sub generateBlasts {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;

  my $nbBlastsToDo;
  my $nbBlastsToDoPerSample = LoadConfig::getParam($rH_cfg, 'generateBlasts', 'blastToDoPerSample');
  my $nbBlastsToDoPerLane = LoadConfig::getParam($rH_cfg, 'generateBlasts', 'blastToDoPerLane');

  if (!defined($nbBlastsToDoPerSample) || ($nbBlastsToDoPerSample eq "") || ($nbBlastsToDoPerSample < 1) || (ref($nbBlastsToDoPerSample) eq "ARRAY" && scalar(@{$nbBlastsToDoPerSample}) == 0)) {
    $nbBlastsToDo = ceil($nbBlastsToDoPerLane / (scalar (@$rAoH_sample)));
  } else {
    $nbBlastsToDo = $nbBlastsToDoPerSample;
  }

  my $dependencies = getDependencies($step, $rH_cfg);

  for my $rH_sample (@$rAoH_sample) {
    my $outputPrefix = $runDirectory.'/' . $UNALIGNED_DIR . '.'.$lane.'/Blast_sample/'.$rH_sample->{'processingSheetId'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'};
    my $read1File = getFastqFilename($runDirectory, $lane, $rH_sample, 1);
    my $read2File = getFastqFilename($runDirectory, $lane, $rH_sample, 2);
    my $ro_job = new Job();
    
    if ($nbReads == 1) {
      $ro_job->testInputOutputs([$read1File],[$outputPrefix.'.R1.RDP.blastHit_20MF_species.txt']);
    } else {
      $ro_job->testInputOutputs([$read1File, $read2File],[$outputPrefix.'.R1.RDP.blastHit_20MF_species.txt']);
    }

    if (!$ro_job->isUp2Date()) {
      my $command = 'mkdir -p '.$runDirectory.'/' . $UNALIGNED_DIR . '.'.$lane.'/Blast_sample && ';
      $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'generateBlasts','moduleVersion.tools') . ' && ' ;

      if($nbReads == 1) {
        $command .= 'runBlast.sh '.$nbBlastsToDo.' '. $outputPrefix.' '. $read1File;
      } else {
        $command .= 'runBlast.sh '.$nbBlastsToDo.' '. $outputPrefix.' '. $read1File . ' ' . $read2File;
      }
      $ro_job->addCommand($command);
      my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "generateBlasts", "$runID.$lane", 'blast_'.$rH_sample->{'processingSheetId'}, $dependencies, $rH_sample->{'processingSheetId'}, $ro_job);
      push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
    }
  }
  return;
}

sub align {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;
  my $rAoH_sampleLanes  = $rAoH_sample;

  my $jobDependency = getDependencies($step, $rH_cfg);

  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $sampleName = $rH_laneInfo->{'name'};
    my $processingId = $rH_laneInfo->{'processingSheetId'};
    $step->{'jobIds'}->{$processingId} =();
    my $libSource = $rH_laneInfo->{'libSource'}; # gDNA, cDNA, ...
    my $ref = getGenomeReference($rH_laneInfo->{'referenceMappingSpecies'}, $rH_laneInfo->{'ref'}, $libSource, 'bwa');
    if (!defined($ref)) {
      print STDERR "Skipping alignment for sample '$sampleName'; No reference genome found for species '". (defined($rH_laneInfo->{'ref'}) ? $rH_laneInfo->{'ref'} : ""). "'.\n";
      next;
    }

    my $pair1 = getFastqFilename($runDirectory, $lane, $rH_laneInfo, 1);
    my $pair2;
    if ($nbReads > 1) {
      $pair2 = getFastqFilename($runDirectory, $lane, $rH_laneInfo, 2);
    }

    my $rgId = $rH_laneInfo->{'libraryBarcode'} . "_" . $runID . "_" . $rH_laneInfo->{'lane'};
    my $rgSampleName = $rH_laneInfo->{'name'};
    my $rgLibrary = $rH_laneInfo->{'libraryBarcode'};
    my $rgPlatformUnit = $runID . "_" . $rH_laneInfo->{'lane'};
    my $rgCenter = LoadConfig::getParam( $rH_cfg, 'mem', 'bwaInstitution' );

    my $outputAlnDir = $runDirectory.'/'. $ALIGNED_DIR . '.'.$lane.'/alignment/'.$sampleName .'/run' .$runID . "_" . $rH_laneInfo->{'lane'};

    my $outputAlnPrefix = $outputAlnDir.'/'.$sampleName.'.'.$rH_laneInfo->{'libraryBarcode'};

    my $ro_bwaJob = BWA::mem($rH_cfg, $sampleName, $pair1, $pair2, $pair1, $outputAlnPrefix, $rgId, $rgSampleName, $rgLibrary, $rgPlatformUnit, $rgCenter, $ref);
    if(!$ro_bwaJob->isUp2Date()) {
      print 'mkdir -p '.$outputAlnDir."\n";
      my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mem", "$runID.$lane", 'BWA_MEM_'.$processingId, $jobDependency, $processingId, $ro_bwaJob);
      push (@{$step->{'jobIds'}->{$processingId}}, $jobId);
    }
  }
  return;
}

sub laneMetrics {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;

  my $rAoH_sampleLanes  = $rAoH_sample;

  my $first=1;
  my %downloadedBedFiles;
  
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $libSource = $rH_laneInfo->{'libSource'}; # gDNA, cDNA, ...
    my $ref = getGenomeReference($rH_laneInfo->{'referenceMappingSpecies'}, $rH_laneInfo->{'ref'}, $libSource, 'fasta');
    if (!defined($ref)) {
      #skipped alignment
      next;
    }

    my $sampleName = $rH_laneInfo->{'name'};
    my $processingId = $rH_laneInfo->{'processingSheetId'};
    my $jobDependency = getDependencies($step, $rH_cfg, $processingId);
    $step->{'jobIds'}->{$processingId} =();

    my $directory = $runDirectory.'/' . $ALIGNED_DIR. '.'.$lane.'/alignment/'.$sampleName."/run".$runID."_".$rH_laneInfo->{'lane'}."/";
    my $sortedLaneBamFile = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.bam';
    my $sortedLaneDupBamFile = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.dup.bam';
    my $outputMetrics = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.dup.metrics';

    my $rO_job = Picard::markDup($rH_cfg, $sampleName, $sortedLaneBamFile, $sortedLaneDupBamFile, $outputMetrics);
    if(!$rO_job->isUp2Date()) {
      my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "markDup",$runID . "." . $rH_laneInfo->{'lane'}, 'LANEMARKDUP_'.$processingId, $jobDependency, $processingId, $rO_job);
      push (@{$step->{'jobIds'}->{$processingId}}, $jobId);
      push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
    }

    $outputMetrics = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.metrics';
    my $rO_collectMetricsJob = Picard::collectMetrics($rH_cfg, $sortedLaneBamFile, $outputMetrics, $ref);
    if(!$rO_collectMetricsJob->isUp2Date()) {
      my $jobId2 = SubmitToCluster::printSubmitCmd($rH_cfg, "collectMetrics", $runID . "." . $rH_laneInfo->{'lane'}, 'COLLECTMETRICS_'.$processingId, $jobDependency, $processingId, $rO_collectMetricsJob);
      push (@{$step->{'jobIds'}->{$processingId}}, $jobId2);
      push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId2);
    }
    
    $outputMetrics = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.metrics.targetCoverage.txt';
    my $coverageBED = defined(BVATools::resolveSampleBED($rH_cfg, $rH_laneInfo)) ? $runDirectory . "/". BVATools::resolveSampleBED($rH_cfg, $rH_laneInfo) : undef;
    # download bed files?
    if (LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'fetchBedFiles')) {
      for my $bedFile (@{$rH_laneInfo->{'bedFiles'}}) {
        if (!defined($downloadedBedFiles{$bedFile})) {
          print formatCommand("config" => $rH_cfg, "command" => LoadConfig::getParam($rH_cfg, 'depthOfCoverage', 'fetchBedFileCommand'), "runDirectory" => $runDirectory, "runID" => $runID, "lane" => $lane , "bedFile" => $bedFile) . "\n";
          print "\n";
          $downloadedBedFiles{$bedFile} = 1;
        }
      }
    }

    my $rO_coverageJob = BVATools::depthOfCoverage($rH_cfg, $sortedLaneBamFile, $outputMetrics, $coverageBED, $ref);
    if(!$rO_coverageJob->isUp2Date()) {
      my $jobId3 = SubmitToCluster::printSubmitCmd($rH_cfg, "depthOfCoverage", $runID . '.' . $rH_laneInfo->{'lane'}, 'LANEDEPTHOFCOVERAGE_'.$processingId, $jobDependency, $processingId, $rO_coverageJob);
      push (@{$step->{'jobIds'}->{$processingId}}, $jobId3);
      push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId3);
    }

    my ($refDict) = $ref =~ /^(.+)\.[^.]*$/;
    $refDict .= '.dict';
    if(defined($coverageBED)) {
      my ($coverageIntList) = $coverageBED =~ /^(.+)\.[^.]+$/;
      $coverageIntList .= '.interval_list';

      my $rO_intListJob = Tools::generateIntervalList($rH_cfg, $refDict, $coverageBED, $coverageIntList);
      if(!$rO_intListJob->isUp2Date() && !defined($generatedIntervalLists{$coverageIntList})) {
        print $rO_intListJob->getCommand()."\n";
        $generatedIntervalLists{$coverageIntList} = 1;
      }

      $outputMetrics = $directory.$rH_laneInfo->{'name'}.'.'.$rH_laneInfo->{'libraryBarcode'}.'.sorted.metrics.onTarget.txt';
      my $rO_hsMetricsJob = Picard::calculateHSMetrics($rH_cfg, $sortedLaneBamFile, $coverageIntList, $outputMetrics, $ref, $refDict);
      if(!$rO_hsMetricsJob->isUp2Date()) {
        my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "calculateHSMetrics", $runID . '.' . $rH_laneInfo->{'lane'}, 'CALCULATEHSMETRICS_'.$processingId, $jobDependency, $processingId, $rO_hsMetricsJob);
        push (@{$step->{'jobIds'}->{$processingId}}, $jobId);
        push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
      }
    }
  }
  return;
}

sub generateBamMd5 {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;

  for my $rH_sample (@$rAoH_sample) {
    my $sampleName = $rH_sample->{'name'};
    my $processingId = $rH_sample->{'processingSheetId'};
    my $libSource = $rH_sample->{'libSource'}; # gDNA, cDNA, ...
    my $ref = getGenomeReference($rH_sample->{'referenceMappingSpecies'}, $rH_sample->{'ref'}, $libSource, 'bwa');
    if (!defined($ref)) {
      #skipped alignment
      next;
    }

    my $dependencies = getDependencies($step, $rH_cfg, $processingId);
    my $directory = $runDirectory.'/' . $ALIGNED_DIR . '.'.$lane.'/alignment/'.$sampleName."/run".$runID."_".$lane."/";
    my $sortedLaneFile = $directory . $sampleName.'.'.$rH_sample->{'libraryBarcode'}.'.sorted';

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$sortedLaneFile.'.bai'],[$sortedLaneFile.'.bai.md5']);
    if (!$ro_job->isUp2Date()) {
      my $command = 'md5sum -b '.$sortedLaneFile . '.bai > ' . $sortedLaneFile . '.bai.md5';
      $ro_job->addCommand($command);
      my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "generateBamMD5", "$runID.$lane", 'bmd5_'.$processingId, $dependencies, $processingId, $ro_job);
      push (@{$step->{'jobIds'}->{$processingId}}, $jobId);
      push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
    }
  }
  return;
}

sub startCopyNotification {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;

  my $dependencies = getDependencies($step, $rH_cfg);

  if (LoadConfig::getParam($rH_cfg, 'startCopyNotification', 'sendNotification')) {
    my $command = formatCommand("config" => $rH_cfg, "command" => LoadConfig::getParam($rH_cfg, 'startCopyNotification', 'notificationCommand'), "runDirectory" => $runDirectory, "runID" => $runID, "lane" => $lane);

    my $ro_job = new Job();
    $ro_job->addCommand($command);

    my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "startCopyNotification", "$runID.$lane", 'startCopyNotification', $dependencies, undef, $ro_job);
    push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
  }
  return;
}

sub copy {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;
  my $rAoH_sample    = shift;

  my $dependencies = getDependencies($step, $rH_cfg);
  my $ro_job = new Job();
  my $destinationFolder;

  $destinationFolder= LoadConfig::getParam($rH_cfg, 'copy','destinationFolder');

  my $excludeFastq = "";
  # Check if sample has a reference genome and therefore associated BAM files
  for my $rH_sample (@$rAoH_sample) {
    my $libSource = $rH_sample->{'libSource'}; # gDNA, cDNA, ...
    my $ref = getGenomeReference($rH_sample->{'referenceMappingSpecies'}, $rH_sample->{'ref'}, $libSource, 'bwa');
    if (defined($ref)) {
      # BAM have been created, therefore no need to rsync fasta files
      $excludeFastq .= " \\\n  --exclude '" . getFastqFilename($runDirectory, $lane, $rH_sample, 1) . "'";
      if ($nbReads > 1) {
        $excludeFastq .= " \\\n  --exclude '" . getFastqFilename($runDirectory, $lane, $rH_sample, 2) . "'";
      }
    }
  }

  my $command = formatCommand("config" => $rH_cfg, "command" => LoadConfig::getParam($rH_cfg, 'copy', 'copyCommand'), "runDirectory" => $runDirectory, "runID" => $runID, "lane" => $lane, "destinationFolder" => $destinationFolder, "excludeFastq" => $excludeFastq);

  $ro_job->addCommand($command);

  my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "copy", "$runID.$lane", 'copy_', $dependencies, undef, $ro_job);
  push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
  return;
}

sub endCopyNotification {
  my $step           = shift;
  my $rH_cfg         = shift;
  my $runDirectory   = shift;
  my $runID          = shift;
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $nbReads        = shift;

  my $dependencies = getDependencies($step, $rH_cfg);

  if (LoadConfig::getParam($rH_cfg, 'endCopyNotification', 'sendNotification')) {
    my $command = formatCommand("config" => $rH_cfg, "command" => LoadConfig::getParam($rH_cfg, 'endCopyNotification', 'notificationCommand'), "runDirectory" => $runDirectory, "runID" => $runID, "lane" => $lane);

    my $ro_job = new Job();
    $ro_job->addCommand($command);

    my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, "endCopyNotification", "$runID.$lane", 'endCopyNotification', $dependencies, undef, $ro_job);
    push (@{$step->{'jobIds'}->{$GLOBAL_DEP_KEY}}, $jobId);
  }
  return;
}

###################
# Utility methods #
###################

# Build the default genome list from the configuration file, used to generate alignment from a species regexp
sub getDefaultGenomes {
  my $genomes = shift;
  my %defaultgenomes;

  for my $genome (split('~', $genomes)) {
    my ($regexp, $species, $build) = split(':', $genome);
    $defaultgenomes{$regexp}{"species"} = $species;
    $defaultgenomes{$regexp}{"build"} = $build;
  }

  return \%defaultgenomes;
}


# Parse on disk the genome list to construct a data structure of the style {species}{build}{program}
sub getGenomeList {
  my $rootDir = shift;

  opendir(ROOT_DIR, $rootDir) or exitWithError("Couldn't open directory ".$rootDir);
  my @speciesDirs =  grep { /^[^\.]/ && -d "$rootDir/$_" } readdir(ROOT_DIR);
  closedir(ROOT_DIR);

  my %genomes;

  for my $speciesDir (@speciesDirs) {
    opendir(SPECIES_DIR, "$rootDir/$speciesDir") or exitWithError("Couldn't open directory $rootDir/$speciesDir");
    my @buildDirs = grep { /^[^\.]/ && -d "$rootDir/$speciesDir/$_" } readdir(SPECIES_DIR);
    closedir(SPECIES_DIR);
    for my $buildDir (@buildDirs) {
      my $fasta = "$rootDir/$speciesDir/$buildDir/fasta/bwa/$buildDir.fasta";
      my $fa = "$rootDir/$speciesDir/$buildDir/fasta/bwa/$buildDir.fa";
      if (-r $fasta) {
        $genomes{$speciesDir}{$buildDir}{"bwa"}=$fasta;
      } elsif (-r $fa) {
        $genomes{$speciesDir}{$buildDir}{"bwa"}=$fa;
      } else {
        #print STDERR "Available Genomes Scan: No BWA reference genome found for the build '$buildDir' of the '$speciesDir' species\n";
      }
      $fasta = "$rootDir/$speciesDir/$buildDir/fasta/$buildDir.fasta";
      $fa = "$rootDir/$speciesDir/$buildDir/fasta/$buildDir.fa";
      if (-r $fasta) {
        $genomes{$speciesDir}{$buildDir}{"fasta"}=$fasta;
      } elsif (-r $fa) {
        $genomes{$speciesDir}{$buildDir}{"fasta"}=$fa;
      } else {
        #print STDERR "Available Genomes Scan: No Fasta reference genome found for the build '$buildDir' of the '$speciesDir' species\n";
      }

    }
  }

  return \%genomes;
}

# Get the corresponding genome reference
# if the sample is gDNA, will use the "reference mapping species", if populated;
# the species will be used otherwise
sub getGenomeReference {
  my $ref           = shift;
  my $species       = shift;
  my $librarySource = shift;
  my $program       = shift;

  my $refpath;

  if ($librarySource eq "gDNA" && defined($ref)) {
    # gDNA alignment with BWA
    my ($refSpecies, $build) = split( ',', $ref );

    if (!defined($refSpecies) || !defined($build)){
      return;
    }

    # Trimming leading/trailing spaces
    $refSpecies =~ s/^\s+|\s+$//g;
    $build =~ s/^\s+|\s+$//g;

    # Replacing spaces with '_'
    $refSpecies =~ s/\s/_/g;
    $build =~ s/\s/_/g;

    $refpath =  $rHoH_genomes->{$refSpecies}->{$build}->{$program}
  }

  if (defined($refpath)) {
    return $refpath;
  } else {
    if (defined($species)) {
      # defaulting to a basic alignement with BWA
      for my $defaultGenomeRegexp (keys %$rHoH_defaultGenomes) {
        if ($species =~ /^\s*$defaultGenomeRegexp\s*$/i) {
          my $refSpecies = $rHoH_defaultGenomes->{$defaultGenomeRegexp}->{"species"};
          my $build = $rHoH_defaultGenomes->{$defaultGenomeRegexp}->{"build"};
          $refpath = $rHoH_genomes->{$refSpecies}->{$build}->{$program}
        }
      }
    }
  }

  return $refpath;

}


sub parseRunInfoFile {
  my $fileName = shift;
  
  my $xml     = new XML::Simple( 'ForceArray' => ['Read'] );
  my $runInfo = $xml->XMLin( $fileName );
  my $nbLanes = $runInfo->{'Run'}->{'FlowcellLayout'}->{'LaneCount'};
  my $nbTemplateReads = 0;
  my @AoH_parsedReads;

  # Get read information from runinfo first
  for my $rh_read (@{$runInfo->{'Run'}->{'Reads'}->{'Read'}}) {
    my $nbCycles = $rh_read->{'NumCycles'};
    my $isIndexed = $rh_read->{'IsIndexedRead'};
    push(@AoH_parsedReads, {nbCycles=>$nbCycles , isIndexed=>$isIndexed, nbLanes=>$nbLanes});
    if($isIndexed ne "Y") {
      $nbTemplateReads++;
    }
  }
  return ($nbTemplateReads, \@AoH_parsedReads);
}


# Return the path to the fastq generated by Casava
sub getFastqFilename {
  my $runDirectory = shift;
  my $lane         = shift;
  my $rH_sample    = shift;
  my $readNb       = shift;
  
  return $runDirectory.'/' . $UNALIGNED_DIR . '.'.$lane.'/Project_nanuq/Sample_'.$rH_sample->{'processingSheetId'}.'/'.$rH_sample->{'processingSheetId'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R' . $readNb . '_001.fastq.gz'
}

sub computeIndexLength {
  my $rAoH_readsInfo = shift;

  my $indexLength = 0;
  for my $rH_readInfo (@{$rAoH_readsInfo}) {
    if($rH_readInfo->{'isIndexed'} eq 'Y') {
      $indexLength += $rH_readInfo->{'nbCycles'};
    }
  }
  return $indexLength;
}


sub formatCommand {
  my %params = @_;
  my $command = $params{'command'};

  $params{'technology'} = LoadConfig::getParam($params{'config'}, '', 'technologyName');

  if (defined($params{'runDirectory'})) {
    $params{'runName'} = basename($params{'runDirectory'});
  }
  
  $params{'Unaligned'} = $UNALIGNED_DIR;
  $params{'Aligned'} = $ALIGNED_DIR;

  my $param;
  while(($param) = ($command =~ /\${(\w+)}/g )) {
    if (defined($params{$param})) {
      $command =~ s/\${$param}/$params{$param}/g;
    }
  }
  return $command;
}

sub getDependencies {
  my $step = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;

  if (!defined($sampleName)) {
    $sampleName = $GLOBAL_DEP_KEY;
  }

  my $dependencies = "";
  # Retrieve the list of step parents
  my @A_stepParents = map {$H_steps{$_}} @{$step->{'parentStep'}};

  # Retrieve the list of lists of step parent job IDs if any
  my @AoA_stepParentJobIds = map {defined $_->{'jobIds'}->{$sampleName} ? $_->{'jobIds'}->{$sampleName} : []} @A_stepParents;
  # Flatten this list
  my @A_stepParentJobIds = map {@$_} @AoA_stepParentJobIds;

  # Concatenate all job IDs with cluster dependency separator
  $dependencies = join (LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep'), map {"\$" . $_} @A_stepParentJobIds);
  return $dependencies;
}




####################################################################
# Methods for calculating index mask (and generating sample sheet) #
####################################################################
sub validateBarcodes {
  my $numberOfMismatch = shift;
  my $rAoH_sample      = shift;

  my @indexes;
  my @collisions;

  my $minAllowedDistance = (2 * $numberOfMismatch) + 1;

  for my $rH_sample (@$rAoH_sample) {
    my $currentIndex = $rH_sample->{'index'};
    $currentIndex =~ s/\-//; # remove dual-index "-"
    
    for my $candidateIndex (@indexes) {
      my $distance = distance($currentIndex, $candidateIndex);
      if ($distance < $minAllowedDistance) {
        push (@collisions, "'$currentIndex' and '$candidateIndex'");
      }
    }
    push(@indexes, $currentIndex);
  }

  if (scalar(@collisions) > 0) {
    exitWithError("Barcode collisions: " . join("; ", @collisions));
  }
  return;
}

sub exitWithError {
  my $message = shift;

  # also print the message to the stdout;
  print "echo '$message';\n exit 1";
  die $message;
}


# One liner returning the hamming distance between two strings of the same length
sub distance {
  my $a = shift;
  my $b = shift;
  return ($a ^ $b) =~ tr/\001-\255//;
}

sub getMask {
  my $lane           = shift;
  my $rAoH_readsInfo = shift;
  my $firstIndex     = shift;
  my $lastIndex      = shift;
  
  my ($rA_headers, $rAoA_sampleSheetDatas) = getSampleSheetContent();
  my $rAoA_laneData = getLaneSampleDatas($lane, $rA_headers, $rAoA_sampleSheetDatas);
  my $rA_laneIdxLengths = getSmallestIndexLength($rAoH_readsInfo, $rA_headers, $rAoA_laneData);

  my $mask = "";
  my $readIndex = 0;
  my $nbTotalIndexBaseUsed = 0;
  
  for my $rH_readInfo (@{$rAoH_readsInfo}) {
    my $nbCycles = $rH_readInfo->{'nbCycles'};
    if (length($mask) != 0) {
      $mask .= ',';
    }

    if ($rH_readInfo->{'isIndexed'} eq "Y") {
      if ($nbCycles >= $rA_laneIdxLengths->[$readIndex]) {
        if ($rA_laneIdxLengths->[$readIndex] == 0 || $lastIndex <= $nbTotalIndexBaseUsed) {
          # Don't use any index base for this read
          $mask .= 'n'.$nbCycles;
        } else {
          my $nbNPrinted = 0;
          my $nbIndexBaseUsed = 0;
          
          # Ns in the beginning of the index read
          if ($firstIndex > ($nbTotalIndexBaseUsed + 1)) {
            $nbNPrinted = min($nbCycles, $firstIndex-$nbTotalIndexBaseUsed-1);
            if ($nbNPrinted >= $rA_laneIdxLengths->[$readIndex]) {
              $nbNPrinted = $nbCycles;
            }
            $mask .= 'n' . $nbNPrinted;
          }
          
          # Calculate the number of index bases
          $nbIndexBaseUsed = max($rA_laneIdxLengths->[$readIndex] - $nbNPrinted, 0);
          $nbIndexBaseUsed = min($lastIndex-$nbTotalIndexBaseUsed-$nbNPrinted, $nbIndexBaseUsed);
          $nbTotalIndexBaseUsed += $nbIndexBaseUsed + min($nbNPrinted, $rA_laneIdxLengths->[$readIndex]);
          if ($nbIndexBaseUsed > 0) {
            $mask .= 'I'.$nbIndexBaseUsed;
          }
          
          # Ns at the end of the index read
          if (($nbCycles-$nbIndexBaseUsed-$nbNPrinted) > 0) {
            $mask .= 'n'.($nbCycles-$nbIndexBaseUsed-$nbNPrinted);
          }
        }
      } else {
        exitWithError("Cycles for index don't match on lane: ".$lane);
      }
      $readIndex++;
    } else {
      $mask .= 'Y'.$nbCycles;
    }
  }
  return $mask;
}

sub min {
  my $a = shift;
  my $b = shift;
  
  return ($a < $b) ? $a : $b;
}

sub max {
  my $a = shift;
  my $b = shift;
  
  return ($a < $b) ? $b : $a;
}

#Returns the sampleSheeet column header as a reference to an array and the sampleSheet sample data line as a reference to an array of array
sub getSampleSheetContent{
  my @dataLines;
  my @headers;

  #get the sample from nanuq
  my $SSHEET;
  open( $SSHEET, '<', $casavaSheet ) or exitWithError("Can't open sample sheet: " . $casavaSheet);

  #parse header line
  my $header=<$SSHEET>;
  @headers=split( /,/ , $header);

  while ( my $line = <$SSHEET> ) {
    #data line
    my @values = split( /,/, $line );
    if(@values != @headers){
      exitWithError("Missing data columns in sampleSheet for row: " . join(", " , @values));
    }
    push(@dataLines, \@values );
  }
  close($SSHEET);
  return (\@headers, \@dataLines);
}


# Returns an array reference of array reference of all the lane data;
sub getLaneSampleDatas{
  my $lane                  = shift;
  my $rA_headers            = shift;
  my $rAoA_sampleSheetDatas = shift;
  my @reply;

  my $laneColumnIdx = getColumnHeaderIndex('Lane', $rA_headers);
  for my $rA_values (@$rAoA_sampleSheetDatas ) {
    my $tempLane = $rA_values->[$laneColumnIdx];
    if($tempLane == $lane){
      push(@reply, $rA_values);
    }
  }
  return \@reply;
}

sub getColumnHeaderIndex{
  my $columnName = shift;
  my $rA_headers = shift;

  for (my $idx = 0 ; $idx < @$rA_headers ; $idx++ ) {
    if($rA_headers->[$idx] eq $columnName){
      return $idx;
    }
  }
  return -1;
}

# Returns an array reference of the smallest index length for the given lane data;
sub getSmallestIndexLength{
  my $rAoH_readsInfo  = shift;
  my $rA_headers      = shift;
  my $rAoA_sampleData = shift;

  my @runIdxLengths;

  for my $rH_readInfo (@{$rAoH_readsInfo}) {
    if($rH_readInfo->{'isIndexed'} eq 'Y') {
      push(@runIdxLengths, $rH_readInfo->{'nbCycles'});
    }
  }

  my $maxSampleIndexRead = 0;
  # No indexes read on sequencer
  if(@runIdxLengths > 0) {
    # find smallest index per index-read per lane
    my $indexColumnIdx = getColumnHeaderIndex('Index', $rA_headers);

    for my $rA_values (@$rAoA_sampleData ) {
      my @libraryIndexes = split('-', $rA_values->[$indexColumnIdx]);

      for(my $idx=0; $idx < @libraryIndexes; $idx++) {
        $maxSampleIndexRead = $idx if ($idx > $maxSampleIndexRead);
        if(length($libraryIndexes[$idx]) > 0) {
          if(length($libraryIndexes[$idx]) < $runIdxLengths[$idx]) {
            $runIdxLengths[$idx] = length($libraryIndexes[$idx]);
          }
        }
      }
    }
  }
  elsif(@$rAoA_sampleData > 1) {
    exitWithError("Multiple samples on a lane, but no indexes were read from the sequencer.");
  }

  # In the case of single-index lane in a dual index run
  for (my $idx=$maxSampleIndexRead+1; $idx < @runIdxLengths; $idx++) {
    $runIdxLengths[$idx] = 0;
  }

  return \@runIdxLengths;
}


#Generates the lane samplesheet on disk and return an array of sampleInfo with the indexToUse initialized.
sub generateIlluminaLaneSampleSheet {
  my $lane           = shift;
  my $runDirectory   = shift;
  my $rAoH_readsInfo = shift;
  my $mask           = shift;
  
  my @readMask = split(',', $mask);
  
  #init samplesheet data
  my ($rA_headers, $rAoA_sampleSheetDatas) = getSampleSheetContent($casavaSheet);
  my $indexColumnIdx = getColumnHeaderIndex('Index', $rA_headers);
  my $laneColumnIdx = getColumnHeaderIndex('Lane', $rA_headers);
  my $sampleIDColumnIdx = getColumnHeaderIndex('SampleID', $rA_headers);
  my $sampleRefIdx = getColumnHeaderIndex('SampleRef',$rA_headers);

  #validate presence of mandatory columns
  if ( $indexColumnIdx==-1  || $laneColumnIdx==-1  || $sampleIDColumnIdx==-1 || $sampleRefIdx==-1 )  {
    exitWithError("Missing header columns");
  }

  my @retVal;
  #print smallest index to use in lane file instead of index found in sampleSheet


  #print sample header
  my $file = $runDirectory.'/nanuqSampleSheet.'.$lane.'.csv';
  my $fh;
  open( $fh, '>', $file ) or exitWithError("Can't write nanuq sample sheet: " . $file);
  print {$fh} join(',', @$rA_headers);

  my $rAoA_laneData = getLaneSampleDatas($lane, $rA_headers, $rAoA_sampleSheetDatas);
  my $rA_laneIndexLength=getSmallestIndexLength($rAoH_readsInfo, $rA_headers, $rAoA_laneData);
  my $laneHasOneSample= @$rAoA_laneData==1?1:0;
  my $sampleAreMixed=areSamplesIndexMixed($rA_headers, $rAoA_laneData);
  #print STDERR "Sample are mixed in Lane: $lane is : " . $sampleAreMixed . " Lane has one sample is: " . $laneHasOneSample ."\n";

  for my $rA_values (@$rAoA_laneData ) {
    my $indexToUse = "";
    #sample barcode
    print {$fh} $rA_values->[0];

    for ( my $idx = 1 ; $idx < @$rA_values ; $idx++ ) {
      print {$fh} ',';
      my $columnValue=$rA_values->[$idx];

      if ( $idx == $indexColumnIdx ) {
        if (length($columnValue) > 0) {
          if (!$laneHasOneSample){
            # index to use for lane with more than one sample
            
            my @sampleIndexes = split('-', $columnValue);
            my $nbIndex=@$rA_laneIndexLength;

            if ($sampleAreMixed){
              # we have a mixed of index in the sample, there are sample with one or 2 index, ignore the second index in the samplesheet
              $nbIndex=1;
            }

            for (my $indexIdx=0; $indexIdx < $nbIndex; $indexIdx++) {
              my $currentReadMask = @readMask[($indexIdx > 0) ? 2 : 1];
              my $numberOfIgnoredLeadingIndexBases = 0;
              my $numberOfIndexBasesToUse=0;
              if ($currentReadMask =~ /(n\d+)?(I\d+)(n\d+)?/) {
                if (defined($1)) {
                  $numberOfIgnoredLeadingIndexBases = substr($1, 1);
                }
                if (defined($2)) {
                  $numberOfIndexBasesToUse = substr($2, 1);
                }
              }
              # trim index to smallest lane index
              if ($indexIdx < @sampleIndexes){
                # the sample has this index
                my $index = substr( $sampleIndexes[$indexIdx], $numberOfIgnoredLeadingIndexBases, $numberOfIndexBasesToUse);
                if ($indexIdx > 0 && length($index) > 0) {
                  $indexToUse .= '-';
                }
                $indexToUse .= $index;
              }
            }
          }
          #print indexToUse in samplesheet
          print {$fh} $indexToUse;
        }
      } else {
        #other sample data
        print {$fh} $columnValue;
      }
    } # for sample column

    my %sampleInfo;
    
    if ($indexToUse eq "") {
      $indexToUse = "NoIndex";
    }
    
    $sampleInfo{'processingSheetId'} = $rA_values->[$sampleIDColumnIdx];
    $sampleInfo{'lane'} = $rA_values->[$laneColumnIdx];
    $sampleInfo{'ref'} = $rA_values->[$sampleRefIdx];
    $sampleInfo{'index'} = $indexToUse;

    push(@retVal, \%sampleInfo);

    #close lane sample sheet

  } # for each samples
  close($fh);
  return \@retVal;
}

# Returns 1 if sample index are mixed otherwise return 0;
sub areSamplesIndexMixed{
  my $rA_headers=shift;
  my $rAoA_sampleData=shift;

  my $indexColumnIdx = getColumnHeaderIndex('Index', $rA_headers);

  my $previousIndexSize=-1;
  for my $rA_values (@$rAoA_sampleData ) {
    my @libraryIndexes = split('-', $rA_values->[$indexColumnIdx]);
    if($previousIndexSize == -1){
      $previousIndexSize=@libraryIndexes;
    } else {
      if( $previousIndexSize != @libraryIndexes){
        return 1;
      }
    }
  }
  return 0;
}

1;
