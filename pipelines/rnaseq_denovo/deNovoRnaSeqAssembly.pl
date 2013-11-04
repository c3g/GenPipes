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

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin/../../lib";

# Dependencies
#-----------------
use Getopt::Std;
use LoadConfig;
#use Data::Dumper;
use SampleSheet;
use Trinity;
use SubmitToCluster;

# Globals
#---------------------------------
my @steps = (
  {
    'name'   => 'normalization',
    'loop'   => 'global',
    'parent' => undef
  },
  {
    'name'   => 'rnaSeqDeNovoAssembly',
    'loop'   => 'global',
    'parent' => 'normalization'
  },
  {
    'name'   => 'abundance',
    'loop'   => 'sample',
    'parent' => 'rnaSeqDeNovoAssembly'
  }
);

my %globalDep;
foreach my $step (@steps) {
  $globalDep{$step->{'name'}} = {};
}

my %groupDone;

main();

# SUB
#-----------------
sub getUsage {
  my $usage = <<END;
Usage: perl $0 -h | -c FILE -s number -e number -n FILE [-w DIR]
  -h  help and usage
  -c  .ini config file
  -s  start step, inclusive
  -e  end step, inclusive
  -n  nanuq sample sheet
  -w  work directory (default current)

Steps:
END

  # List and number step names
  for (my $i = 1; $i <= @steps; $i++) {
    $usage .= $i . "- " . $steps[$i - 1]->{'name'} . "\n";
  }

  return $usage;
}

sub main {
  # Check options
  my %opts;
  getopts('hc:s:e:n:w:', \%opts);

  if (defined($opts{'h'}) ||
     !defined($opts{'c'}) ||
     !defined($opts{'s'}) ||
     !defined($opts{'e'}) ||
     !defined($opts{'n'})) {
    die (getUsage());
  }

  my $configFile = $opts{'c'};
  my $startStep = $opts{'s'};
  my $endStep = $opts{'e'};
  my $nanuqSampleSheet = $opts{'n'};
  my $workDirectory = $opts{'w'};

  # Get config and sample values
  my %cfg = LoadConfig->readConfigFile($configFile);
  my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($nanuqSampleSheet);

  SubmitToCluster::initPipeline($workDirectory);

  foreach my $step (@steps) {
    my $stepName = $step->{'name'};
    my $rSub_step = \&$stepName;
    my $stepLoop = $step->{'loop'};

    if ($stepLoop eq 'sample') {
      
    } elsif ($stepLoop eq 'global') {
      my $jobIdVar = &$rSub_step($step, \%cfg, $rHoAoH_sampleInfo);
      $globalDep{$rSub_step}->{'global'} = $jobIdVar;
    }
  }
}

sub normalization {
  my $depends = shift;
  my $rH_cfg = shift;

  my $rO_job = Trinity::normalize_by_kmer_coverage($rH_cfg);
  SubmitToCluster::printSubmitCmd($rH_cfg, "normalization", undef, 'NORMALIZATION', undef, undef, $rO_job);
}

sub rnaSeqDeNovoAssembly {
  my $depends = shift;
  my $rH_cfg = shift;

  my $rO_job = Trinity::trinity($rH_cfg);
  SubmitToCluster::printSubmitCmd($rH_cfg, "deNovoAssembly_trinity_no_butterfly", undef, 'TRINITY_NO_BUTTERFLY', undef, undef, $rO_job);

}

sub abundance {
  my $depends = shift;
  my $rH_cfg = shift;

  my $rO_job = Trinity::abundance($rH_cfg, "$workDirectory/results/assembly/Trinity.fasta", "rsem", $leftList, $rightList);
  SubmitToCluster::printSubmitCmd($rH_cfg, "abundance_rsem", undef, 'ABUNDANCE_RSEM', undef, undef, $rO_job);
}


sub deNovoAssembly_old {
  my $depends            = shift;
  my $rH_cfg             = shift;
  my $sampleName         = shift;
  my $rAoH_sampleLanes   = shift;
  my $rHoH_groupInfo     = shift;
  my $rH_aliasSampleInfo = shift;
  my $jobDependency      = undef;
  my $group              = $rH_aliasSampleInfo->{$sampleName}{'group_name'};
  if ($depends > 0) {
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
