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

B<Joel Fillon> - I<joel.fillon@mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing


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

ParallelBlast.pl

Parallelize

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
#-----------------
use Getopt::Std;
use LoadConfig;
use SampleSheet;
use SubmitToCluster;
use Trinity;

# Globals
#---------------------------------
my @steps = (
  {
    'name'   => 'normalization',
    'loop'   => 'global',
    'parent' => undef
  },
  {
    'name'   => 'trinity',
    'loop'   => 'global',
    'parent' => 'normalization'
  },
  {
    'name'   => 'trinityQC',
    'loop'   => 'global',
    'parent' => 'trinity'
  },
  {
    'name'   => 'blastSplitQuery',
    'loop'   => 'global',
    'parent' => 'trinity'
  },
  {
    'name'   => 'blast',
    'loop'   => 'global',
    'parent' => 'blastSplitQuery'
  },
  {
    'name'   => 'blastMergeResults',
    'loop'   => 'global',
    'parent' => 'blast'
  },
  {
    'name'   => 'rsemPrepareReference',
    'loop'   => 'global',
    'parent' => 'trinity'
  },
  {
    'name'   => 'rsem',
    'loop'   => 'sample',
    'parent' => 'rsemPrepareReference'
  },
  {
    'name'   => 'DGE',
    'loop'   => 'global',
    'parent' => 'rsem'
  }
);


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

  for (my $i = $startStep; $i <= $endStep; $i++) {
    my $step = $steps[$i - 1];
    my $stepName = $step->{'name'};
    $step->{'jobIds'} = ();

    if ($step->{'loop'} eq 'sample') {
      foreach my $sample (keys %$rHoAoH_sampleInfo) {
        &$stepName(\%cfg, $step, $workDirectory, $sample);
      }
    } else {  # Global step
      &$stepName(\%cfg, $step, $workDirectory);
    }
  }
}

sub getStepParent {
  my $stepParentName = shift;

  my $stepParent = undef;
  if (defined($stepParentName)) {
    foreach my $step (@steps) {
      if ($step->{'name'} eq $stepParentName) {
        $stepParent = $step;
      }
    }
  }
  return $stepParent;
}

sub submitJob {
  my $rH_cfg = shift;
  my $step = shift;
  my $sample = shift;
  my $rO_job = shift;

  my $stepName = $step->{'name'};
  my $jobIdPrefix = uc($stepName);
  my $jobIds = $jobIdPrefix . "_JOB_IDS";
  if (defined $sample) {
    $jobIdPrefix .= "_" . $sample;
  }

  my $dependencies = "";

  my $stepParentName = $step->{'parent'};
  my $stepParent = getStepParent($stepParentName);
  if (defined($stepParent)) {
    $dependencies = join (":", map {"\$" . $_} @{$stepParent->{'jobIds'}});
  }

  my $jobId = SubmitToCluster::printSubmitCmd($rH_cfg, $stepName, undef, $jobIdPrefix, $dependencies, $sample, $rO_job);

  push (@{$step->{'jobIds'}}, $jobId);
}

# Return parameter value from configuration file or die if this parameter is not defined
sub getParam {
  my $rH_cfg = shift;
  my $section = shift;
  my $paramName = shift;

  my $paramValue = LoadConfig::getParam($rH_cfg, $section, $paramName);

  if ($paramValue) {
    return $paramValue;
  } else {
    die "Error: parameter \"[" . $section . "] " . $paramName . "\" is not defined in the configuration .ini file";
  }
}

# Return module load command string from the given list of modules as [section, moduleVersion] from configuration file
sub moduleLoad {
  my $rH_cfg = shift;
  my $rA_modules = shift;

  return "module load " . join(" ", map {getParam($rH_cfg, $_->[0], $_->[1])} @$rA_modules) . "\n";
}


# Step sub
#---------

sub normalization {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $rO_job = Trinity::normalize_by_kmer_coverage($rH_cfg, $workDirectory);
  submitJob($rH_cfg, $step, undef, $rO_job);
}

sub trinity {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $rO_job = Trinity::trinity($rH_cfg, $workDirectory);
  submitJob($rH_cfg, $step, undef, $rO_job);
}

sub trinityQC {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $rO_job = Trinity::trinityQC($rH_cfg, $workDirectory);
  submitJob($rH_cfg, $step, undef, $rO_job);
}

sub blastSplitQuery {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();
  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    $command .= moduleLoad($rH_cfg, [
      ['blast', 'moduleVersion.exonerate']
    ]);

    my $chunkDir = "\$WORK_DIR/blast/chunks";
    $command .= "mkdir -p $chunkDir\n";
    $command .= "fastasplit -f \$WORK_DIR/trinity_out_dir/Trinity.fasta -o $chunkDir -c " . getParam($rH_cfg, 'blast', 'blastJobs') . " \\\n";

    $rO_job->addCommand($command);
  }
  submitJob($rH_cfg, $step, undef, $rO_job);
}

sub blast {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $numJobs = getParam($rH_cfg, 'blast', 'blastJobs');

  for (my $jobIndex = 0; $jobIndex < $numJobs; $jobIndex++) {
    my $chunkIndex = sprintf("%07d", $jobIndex);

    my $rO_job = new Job();
    if (!$rO_job->isUp2Date()) {
      my $command = "\n";

      $command .= moduleLoad($rH_cfg, [
        ['blast', 'moduleVersion.tools'],
        ['blast', 'moduleVersion.exonerate'],
        ['blast', 'moduleVersion.blast']
      ]);

      my $cores = getParam($rH_cfg, 'blast', 'blastCPUperJob');
      my $program = getParam($rH_cfg, 'blast', 'blastProgram');
      my $db = getParam($rH_cfg, 'blast', 'blastDb');
      my $options = getParam($rH_cfg, 'blast', 'blastOptions');
      my $chunkDir = "\$WORK_DIR/blast/chunks";
      my $chunkQuery = "$chunkDir/Trinity.fasta_chunk_$chunkIndex";
      my $chunkResult = "$chunkDir/$program" . "_Trinity_$db" . "_chunk_$chunkIndex.tsv";

      $command .= "parallelBlast.pl -file $chunkQuery --OUT $chunkResult -n $cores --BLAST \'$program -db $db $options\'";

      $rO_job->addCommand($command);
    }
    submitJob($rH_cfg, $step, "blast_chunk_$jobIndex", $rO_job);
  }
}

sub blastMergeResults {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $numJobs = getParam($rH_cfg, 'blast', 'blastJobs');

  my $rO_job = new Job();
  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    my $program = getParam($rH_cfg, 'blast', 'blastProgram');
    my $db = getParam($rH_cfg, 'blast', 'blastDb');
    my $blastDir = "\$WORK_DIR/blast";
    my $chunkResults = "$blastDir/chunks/$program" . "_Trinity_$db" . "_chunk_*.tsv";
    my $result = "$blastDir/$program" . "_Trinity_$db.tsv";

    $command .= "cat $chunkResults > $result";

    $rO_job->addCommand($command);
  }
  submitJob($rH_cfg, $step, undef, $rO_job);
}

sub rsemPrepareReference {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $rO_job = Trinity::rsemPrepareReference($rH_cfg, $workDirectory);
  submitJob($rH_cfg, $step, undef, $rO_job);
}

sub rsem {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;
  my $sample = shift;

  my $rO_job = Trinity::rsem($rH_cfg, $workDirectory, $sample);
  submitJob($rH_cfg, $step, $sample, $rO_job);
}

sub edgeR {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $rO_job = Trinity::edgeR($rH_cfg, $workDirectory);
  submitJob($rH_cfg, $step, undef, $rO_job);
}
