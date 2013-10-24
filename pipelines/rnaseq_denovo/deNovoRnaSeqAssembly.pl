#!/usr/bin/perl

=head1 NAME

I<deNovoTranscriptomeAssemly>

=head1 SYNOPSIS

deNovoTranscriptomeAssemly.pl B<args> [-f -c -n -s -e]

=head1 DESCRIPTION

B<deNovoTranscriptomeAssemly> Is the main de novo RNA assembly pipeline.

=head1 AUTHORS

B<David Morais> - I<dmorais@cs.bris.ac.uk>

B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing



B<The pipeline specific libs should be in a dir called lib/ placed in the same dir this as this script.>

I<Pipeline specific Libs:>

B<GetFastaAlias>

B<MergeFastq>

B<LoadConfig>

B<SampleSheet>

B<BAMtools>

B<SplitFile>

B<Trinity>

B<BLAST>

B<SampleSheet>

B<SequenceDictionaryParser>

B<SubmitToCluster>

B<Trimmomatic>

B<BWA>

B<HtseqCount>

B<DiffExpression>


=head1 Scripts

This pipeline uses a set of scripts. 
You should create a dir called B<script>
and place it in the same dir as the 
pipeline.

B<This is the list of scripts>

AliasFastqMerge.py

create_gtf.sh

fastalength

fastasplit

filter-duplicates-1.0-jar-with-dependencies.jar

generate_BLAST_HQ.sh

getStat.sh

HtSeq_full_matrix.sh

HtSeq_temp_matrix.sh

ParallelBlast.pl

Parallelize

=cut

# Strict Pragmas
#---------------------
use strict;
use warnings;

#---------------------

BEGIN {
  #Makesure we can find the GetConfig::LoadModules module relative to this script install
  use File::Basename;
  use Cwd 'abs_path';

  my ( undef, $mod_path, undef ) = fileparse( abs_path(__FILE__) );
  unshift @INC, $mod_path . "lib";

}

# Dependencies
#-----------------
use Getopt::Std;
use File::Basename;
use GetFastaAlias;
#use MergeFastq;
use LoadConfig;
use Data::Dumper;
use SampleSheet;
use BAMtools;
use SplitFile;
use Trinity;
use BLAST;
use SampleSheet;
use SequenceDictionaryParser;
use SubmitToCluster;
use BWA;

# Globals
#---------------------------------
my @steps;
push(@steps, {'name' => 'normalization'});
push(@steps, {'name' => 'deNovoAssembly'});

my %groupDone;
my $workDir;

&main();

# SUB
#-----------------
sub printUsage {
  print "\nUsage: perl " . $0 . " args [-f -c -s -e -n]\n";
  print "\t-h  help and usage\n";
  print "\t-c  config file\n";
  print "\t-s  start step, inclusive\n";
  print "\t-e  end step, inclusive\n";
  print "\t-n  nanuq sample sheet\n";
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
  getopts('hc:s:e:n:w:', \%opts);

  if (defined($opts{'h'}) || !defined($opts{'c'}) || !defined($opts{'s'}) || !defined($opts{'e'}) || !defined($opts{'n'}) || !defined($opts{'w'})) {
    printUsage();
    exit(1);
  }

  my %cfg = LoadConfig->readConfigFile($opts{'c'});
  my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($opts{'n'});
  $workDir = abs_path($opts{'w'});

  SubmitToCluster::initPipeline($workDir);

  for (my $currentStep = $opts{'s'} - 1; $currentStep <= ($opts{'e'} - 1); $currentStep++) {
    my $fname = $steps[$currentStep]->{'name'};
    my $subref = \&$fname;

    my $jobIdVar = &$subref($currentStep, \%cfg, $rHoAoH_sampleInfo);
  }
}

sub normalization {
  my $depends = shift;
  my $rH_cfg = shift;

  my $leftList = "$workDir/reads/left_pair1.fastq.gz.list";
  my $rightList = "$workDir/reads/right_pair2.fastq.gz.list";

  print "find $workDir/reads/ -name *pair1.*.fastq.gz > $leftList\n";
  print "find $workDir/reads/ -name *pair2.*.fastq.gz > $rightList\n";

  my $rO_job = Trinity::normalize_by_kmer_coverage($rH_cfg, $workDir, $leftList, $rightList);

  SubmitToCluster::printSubmitCmd($rH_cfg, "normalize", undef, 'NORMALIZE', undef, undef, $rO_job);
}

sub deNovoAssembly {
  my $depends          = shift;
  my $rH_cfg           = shift;
  my $sampleName       = shift;
  my $rAoH_sampleLanes = shift;
  my $rHoH_groupInfo   = shift;
  my $rH_aliasSampleInfo = shift;
  my $jobDependency    = undef;
  my $group            = $rH_aliasSampleInfo->{$sampleName}{'group_name'};
  if ( $depends > 0 ) {
      $jobDependency = '$GROUP_JOB_IDS';
  }
  
   
  # Chrysalis
  #---------------
  for my $rH_laneInfo (@$rAoH_sampleLanes) {

      #next if the one sample of the group is already done
      next if ( exists $groupDone{$group} );
      print "if [ ! -d " . $group . "/output_jobs ] ; then mkdir -p " . $group . "/output_jobs ; fi\n";
      print "if [ ! -d assembly/" . $group . "/output_jobs ]; then  mkdir -p assembly/" . $group . " ; fi\n";
      
      print "CHRYSALIS_JOB_IDS=\"\"\n";

      my $rH_chrysalisDetails = Trinity::chrysalis( $rH_cfg, $group, $rH_laneInfo, $rHoH_groupInfo->{$group}->{'left'} , $rHoH_groupInfo->{$group}->{'right'}); 
      $groupDone{$group} = 1;
      my $chrysalisJobId = undef;

      if ( length( $rH_chrysalisDetails->{'command'} ) > 0 ) {
          $chrysalisJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "trinity", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'CHRYSALIS', $jobDependency, $group, $rH_chrysalisDetails->{'command'} );
          $chrysalisJobId = '$' . $chrysalisJobId;
          print 'CHRYSALIS_JOB_IDS=${CHRYSALIS_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $chrysalisJobId . "\n\n";
          print "GROUP_JOB_IDS=\"\"\n";
      }

      # Split Chrysalis file
      #---------------------------
      print "SPLITCHRYS_JOB_IDS=\"\"\n";
      my $rH_splitChrysalisDetails = SplitFile::splitButterfly( $rH_cfg, $group, $rH_laneInfo );
      my $splitChrysalisJobId = undef;
      if ( length( $rH_splitChrysalisDetails->{'command'} ) > 0 ) {
          $jobDependency = '$CHRYSALIS_JOB_IDS';
          $splitChrysalisJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "splitchrysalis", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'SPLITCHRYSALIS', $jobDependency, $group, $rH_splitChrysalisDetails->{'command'} );
          $splitChrysalisJobId = '$' . $splitChrysalisJobId;
          print 'SPLITCHRYS_JOB_IDS=${SPLITCHRYS_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $splitChrysalisJobId . "\n\n";
      }

      #ButterFly
      #----------
      # This puts the files names in the @Files array (with a 4 digits padding)
      my @files;
      my $nbChunks = LoadConfig::getParam( $rH_cfg, 'butterfly', 'chunks');
      for ( my $i = 0 ; $i < $nbChunks ; $i++ ) {
          push( @files, 'chunk.' . sprintf( "%04d", $i ) . '.txt') ;
      }

      print "BUTTERFLY_JOB_IDS=\"\"\n";

      foreach my $file (@files) {
          my $rH_butterflyDetails = Trinity::butterfly( $rH_cfg, $group, $rH_laneInfo, $file );
          my $butterflyJobId = undef;
          if ( length( $rH_splitChrysalisDetails->{'command'} ) > 0 ) {
              $jobDependency = '$SPLITCHRYS_JOB_IDS';
              $butterflyJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "butterfly", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . '_'. basename($file), 'BUTTERFLY', $jobDependency, $group, $rH_butterflyDetails->{'command'} );
              $butterflyJobId = '$' . $butterflyJobId;
              print 'BUTTERFLY_JOB_IDS=${BUTTERFLY_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $butterflyJobId . "\n\n";

          }

      }

      # Concatenate files and create gft
      #-----------------------------------
      print "CONCAT_JOB_IDS=\"\"\n";
      my $rH_concatDetails = Trinity::concatFastaCreateGtf( $rH_cfg, $group, $rH_laneInfo );
      my $concatJobId = undef;
      if ( length( $rH_splitChrysalisDetails->{'command'} ) > 0 ) {
          $jobDependency = '$BUTTERFLY_JOB_IDS';
          $concatJobId = SubmitToCluster::printSubmitCmd( $rH_cfg, "concat", $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'}, 'CONCAT', $jobDependency, $group, $rH_concatDetails->{'command'} );
          $concatJobId = '$' . $concatJobId;
          print 'CONCAT_JOB_IDS=${CONCAT_JOB_IDS}' . LoadConfig::getParam( $rH_cfg, 'default', 'clusterDependencySep' ) . $concatJobId . "\n\n";

      }

  }

}
