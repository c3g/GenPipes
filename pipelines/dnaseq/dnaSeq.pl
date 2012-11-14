#!/usr/env/perl

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

# Dependencies
#--------------------
use lib '../../lib';
use Getopt::Std;
use LoadConfig;
use LoadModules;
use SampleSheet qw(parseSampleSheet);
use SubmitToCluster;
use BWA;
use MergeBAMs;
#--------------------


# SUB
#--------------------
&main();

my @steps;
push(@steps, {'name' => 'trimAndAlign'});
push(@steps, {'name' => 'mergeLanes'});
push(@steps, {'name' => 'indelRealigner'});
push(@steps, {'name' => 'mergeRealigned'});
push(@steps, {'name' => 'fixmate'});
push(@steps, {'name' => 'sortFixed'});
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

sub printUsage {
  print "\nUsage: perl ".$0." project.csv first_step last_step\n";
  print "\t-c  config file\n";
  print "\t-s  start step, inclusive\n";
  print "\t-e  end step, inclusive\n";
  print "\t-n  nanuq sample sheet\n";
  print "\n";
  print "Steps:\n";
  for(my $idx=1; $idx <= @steps; $idx++) {
    print $idx."- ".$steps[$idx]->{'name'}."\n";
  }
  print "Steps:\n";
}

sub main {
  my %opts = (c=>undef, m=>undef);
  getopts('c:m:s:en:', \%opts);
  
  if (!defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'e'}) || !defined($opts{'n'})) {
    printUsage();
    exit(1);
  }
  
  my %cfg = readConfigFile($opts{'c'});
  my $rHoAoH_sampleInfo = parseSampleSheetAsHash($opts{'n'});

  my $latestBam;
  for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
    my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};

    SubmitToCluster::initSubmit(\%cfg, $sampleName);
    for(my $current = $opts{'s'}; $current <= $opts{'e'}; $current++) {
       my $fname = $steps[$current]->{'name'};
       my $subref = \&$fname;

       # Tests for the first step in the list. Used for dependencies.
       &$subref($current == $opts{'s'}, \%cfg, $sampleName, $rAoH_sampleLanes); 
    }
  }  
}

sub trimAndAlign {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;

  print "BWA_JOB_IDS=\"\"\n";
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $rH_trimDetails = Trimmomatic::trim($rH_cfg, $sampleName, $rH_laneInfo);
    my $trimJobIdVarName=undef;
    if(length($rH_trimDetails->{'command'}) > 0) {
      $trimJobIdVarName = SubmitToCluster::printSubmitCmd($rH_cfg, "trim", undef, 'TRIM', undef, $sampleName, $rH_trimDetails->{'command'});
      #TODO calcReadCounts.sh
    }

    my $rA_commands = BWA::aln($rH_cfg, $sampleName, $rH_laneInfo, $rH_trimDetails->{'pair1'}, $rH_trimDetails->{'pair2'}, $rH_trimDetails->{'single1'}, $rH_trimDetails->{'single2'});
    if(@{$rA_commands} == 3) {
      my $read1JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", undef, 'READ1ALN', $trimJobIdVarName, $sampleName, $rA_commands->{0});
      my $read2JobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", undef, 'READ2ALN', $trimJobIdVarName, $sampleName, $rA_commands->{1});
      my $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", undef, 'BWA', $read1JobId.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$read2JobId, $sampleName, $rA_commands->{2});
      print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$bwaJobId."\n";
    }
    else {
      my $readJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", undef, 'READALN', $trimJobIdVarName, $sampleName, $rA_commands->{0});
      my $bwaJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "aln", undef, 'BWA',  $readJobId, $sampleName, $rA_commands->{1});
      print 'BWA_JOB_IDS=${BWA_JOB_IDS}'.LoadConfig::getParam($rH_cfg, 'aln', 'clusterDependencySep').$bwaJobId."\n";
    } 
  }
}

sub mergeLanes {
  my $depends = shift;
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rAoH_sampleLanes  = shift;
  
  
  my $jobDependency = "";
  if($depends > 0) {
    $jobDependency = '$BWA_JOB_IDS';
  }

  my $command = MergeBAMs::merge($rH_cfg, $sampleName, $rAoH_sampleLanes);
  my $mergeJobId = SubmitToCluster::printSubmitCmd($rH_cfg, "merge", undef, 'MERGELANES', $jobDependency, $sampleName, $command);
}
#push(@steps, {'name' => 'indelRealigner'});
#push(@steps, {'name' => 'mergeRealigned'});
#push(@steps, {'name' => 'fixmate'});
#push(@steps, {'name' => 'sortFixed'});
#push(@steps, {'name' => 'markDup'});
#push(@steps, {'name' => 'metrics'});
#push(@steps, {'name' => 'crestSClip'});
#push(@steps, {'name' => 'sortQname'});
#push(@steps, {'name' => 'countTelomere'});
#push(@steps, {'name' => 'fullPileup'});
#push(@steps, {'name' => 'countTelomere'});
#  print "Step 12: snp and indel calling\n";

1;
