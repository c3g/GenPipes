#!/usr/bin/perl

=head1 NAME

I<dnaSeq>

=head1 SYNOPSIS

dnaSeq.pl

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
push(@steps, {'name' => 'trimAndAlign'});
push(@steps, {'name' => 'mergeLanes'});
push(@steps, {'name' => 'indelRealigner'});
push(@steps, {'name' => 'mergeRealigned'});
push(@steps, {'name' => 'fixmate'});
push(@steps, {'name' => 'markDup'});
push(@steps, {'name' => 'metrics'});
push(@steps, {'name' => 'crestSClip'});
push(@steps, {'name' => 'sortQname'});
push(@steps, {'name' => 'countTelomere'});
push(@steps, {'name' => 'fullPileup'});
push(@steps, {'name' => 'countTelomere'});
#  print "Step 12: snp and indel calling\n";
#  print "Step 13: merge snp calls\n";
#  print "Step 14: filter N streches\n";
#  print "Step 15: flag mappability\n";
#  print "Step 16: snp annotation\n";
#  print "Step 17: snp effect prediction\n";
#  print "Step 18: gene descriptions and GO terms\n";
#  print "Strp 19: dbNSFP annotations\n";
#  print "Step 20: Cosmic annotations\n";


&main();

sub printUsage {
  print "\nUsage: perl ".$0." project.csv first_step last_step\n";
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
  my $rAoH_seqDictionary = SequenceDictionaryParser::readDictFile(\%cfg);

  my $latestBam;
  for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};

    SubmitToCluster::initSubmit(\%cfg, $sampleName);
    for(my $current = $opts{'s'}-1; $current <= ($opts{'e'}-1); $current++) {
       my $fname = $steps[$current]->{'name'};
       my $subref = \&$fname;

       # Tests for the first step in the list. Used for dependencies.
       &$subref($current != ($opts{'s'}-1), \%cfg, $sampleName, $rAoH_sampleLanes, $rAoH_seqDictionary); 
    }
  }  
}

sub trimAndAlign {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  print "BWA_JOB_IDS=\"\"\n";
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo);
    my $trimJobIdVarName=undef;
    if(length($rH_trimDetails->{'command'}) > 0) {
      $trimJobIdVarName = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'TRIM', undef, $sampleName, $rH_trimDetails->{'command'});
      $trimJobIdVarName = '$'.$trimJobIdVarName;
      #TODO calcReadCounts.sh
    }

    my $rA_commands = BWA::aln($rH_cfg, $sampleName, $rH_laneInfo, $rH_trimDetails->{'pair1'}, $rH_trimDetails->{'pair2'}, $rH_trimDetails->{'single1'}, $rH_trimDetails->{'single2'});
    if(@{$rA_commands} == 3) {
      my $read1JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read1.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ1ALN', $trimJobIdVarName, $sampleName, $rA_commands->[0]);
      $read1JobId = '$'.$read1JobId;
      my $read2JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'read2.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READ2ALN', $trimJobIdVarName, $sampleName, $rA_commands->[1]);
      $read2JobId = '$'.$read2JobId;
      my $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'sampe.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA', $read1JobId.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$read2JobId, $sampleName, $rA_commands->[2]);
      $bwaJobId = '$'.$bwaJobId;
      print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$bwaJobId."\n";
    }
    else {
      my $readJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'READALN', $trimJobIdVarName, $sampleName, $rA_commands->[0]);
      $readJobId = '$'.$readJobId;
      my $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", 'samse.'.$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'BWA',  $readJobId, $sampleName, $rA_commands->[1]);
      $bwaJobId = '$'.$bwaJobId;
      print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$bwaJobId."\n";
    } 
  }
  return '$BWA_JOB_IDS';
}

sub mergeLanes {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '$BWA_JOB_IDS';
  }

  my $command = Picard::merge($rH_cfg, $sampleName, $rAoH_sampleLanes);
  my $mergeJobId = undef;
  if(defined($command) && length($command) > 0) {
    $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "merge", undef, 'MERGELANES', $jobDependency, $sampleName, $command);
  }
  return $mergeJobId;
}

sub indelRealigner {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '$MERGELANES_JOB_ID';
  }

  print "mkdir -p $sampleName/realign\n";
  print "REALIGN_JOB_IDS=\"\"\n";
  my $processUnmapped = 1;
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    my $seqName = $rH_seqInfo->{'name'};
    my $command = GATK::realign($rH_cfg, $sampleName, $seqName, $processUnmapped);
    my $intervalJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "indelRealigner", $seqName, 'REALIGN', $jobDependency, $sampleName, $command);
    $intervalJobId = '$'.$intervalJobId;
    print 'REALIGN_JOB_IDS=${REALIGN_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'default', 'clusterDependencySep').$intervalJobId."\n";
    if($processUnmapped == 1) {
      $processUnmapped = 0;
    }
  }
  
  return '${REALIGN_JOB_IDS}';
}

sub mergeRealigned {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${REALIGN_JOB_IDS}';
  }

  my @seqNames;
  for my $rH_seqInfo (@$rAoH_seqDictionary) {
    push(@seqNames, $rH_seqInfo->{'name'});
  }
  my $command = Picard::mergeRealigned($rH_cfg, $sampleName, \@seqNames);
  my $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "mergeRealigned", undef, 'MERGEREALIGN', $jobDependency, $sampleName, $command);
  return $mergeJobId;
}

sub fixmate {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MERGEREALIGN_JOB_ID}';
  }

  my $command = Picard::fixmate($rH_cfg, $sampleName);
  my $fixmateJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "fixmate", undef, 'FIXMATE', $jobDependency, $sampleName, $command);
  return $fixmateJobId;
}

sub markDup {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${FIXMATE_JOB_ID}';
  }

  my $command = Picard::markDup($rH_cfg, $sampleName);
  my $markDupJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "markDup", undef, 'MARKDUP', $jobDependency, $sampleName, $command);
  return $markDupJobId;
}

sub metrics {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  my $rAoH_seqDictionary = shift;

  my $jobDependency = undef;
  if($depends > 0) {
    $jobDependency = '${MARKDUP_JOB_ID}';
  }

  my $command;

  # Collect metrics
  $command = Picard::collectMetrics($rH_cfg, $sampleName);
  my $collectMetricsJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "collectMetrics", undef, 'COLLECTMETRICS', $jobDependency, $sampleName, $command);
  
  # Compute genome coverage
  $command = GATK::genomeCoverage($rH_cfg, $sampleName);
  my $genomeCoverageJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "genomeCoverage", undef, 'GENOMECOVERAGE', $jobDependency, $sampleName, $command);

  # Compute CCDS coverage
#  unless (-e "$sampleName/$sampleName.sorted.dup.CCDS.coverage.done") {
#   my $CCDS_OPTIONS="-T DepthOfCoverage -R $GENOME_FASTA_PATH -I $sampleName/$sampleName.sorted.dup.bam --omitDepthOutputAtEachBase --logging_level ERROR -geneList $REFGENE_FILE --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 100  --start 1 --stop 500 --nBins 499 -dt NONE -L $CCDS_FILE -o $sampleName/$sampleName.sorted.dup.CCDS.coverage --read_buffer_size 17500000";
#   my $CCDS_MOAB_COMMAND="msub -d $CURRENT_DIR -V -l walltime=168:00:0 -q sw -l nodes=1:ppn=6 -j oe -o $output_jobs -N covEx$sampleName -m ae $dependency -m ae -M $EMAIL_NOTIFICATION";
#   print "echo \"rm -f $sampleName/$sampleName.sorted.dup.CCDS.coverage.done ; $JAVA_BIN -Djava.io.tmpdir=$TMP_DIR -Xmx12G -jar $GATK_JAR $CCDS_OPTIONS && touch $sampleName/$sampleName.sorted.dup.CCDS.coverage.done\" | $CCDS_MOAB_COMMAND\n";
#  }

  # Generate IGV track
#  unless (-e "$sampleName/$sampleName.sorted.dup.tdf.done") {
#   my $IGV_MOAB_COMMAND="msub -d $CURRENT_DIR -V -l walltime=48:00:0 -q sw -l nodes=1:ppn=6 -j oe -o $output_jobs -N igv.$sampleName $dependency -m ae -M $EMAIL_NOTIFICATION";
#   print "echo \"rm -f $sampleName/$sampleName.sorted.dup.tdf.done ; $JAVA_BIN -Djava.io.tmpdir=$TMP_DIR -Xmx12G -Djava.awt.headless=true -jar $IGV_TOOLS_JAR count -f min,max,mean $sampleName/$sampleName.sorted.dup.bam $sampleName/$sampleName.sorted.dup.tdf b37 && touch $sampleName/$sampleName.sorted.dup.tdf.done\" | $IGV_MOAB_COMMAND\n";
#  }

  # Compute flags
#  unless (-e "$sampleName/$sampleName.sorted.dup.flagstat.done") {
#   my $FLAGS_MOAB_COMMAND="msub -d $CURRENT_DIR -V -l walltime=48:00:0 -q sw -l nodes=1:ppn=1 -j oe -o $output_jobs -N flag.$sampleName $dependency -m ae -M $EMAIL_NOTIFICATION";
#   print "echo \"rm -f $sampleName/$sampleName.sorted.dup.flagstat.done ; $SAMTOOLS_HOME/samtools flagstat $sampleName/$sampleName.sorted.dup.bam > $sampleName/$sampleName.sorted.dup.flagstat && touch $sampleName/$sampleName.sorted.dup.flagstat.done\" | $FLAGS_MOAB_COMMAND\n";
#  }  
}

#push(@steps, {'name' => 'crestSClip'});
#push(@steps, {'name' => 'sortQname'});
#push(@steps, {'name' => 'countTelomere'});
#push(@steps, {'name' => 'fullPileup'});
#push(@steps, {'name' => 'countTelomere'});
#  print "Step 12: snp and indel calling\n";

1;
