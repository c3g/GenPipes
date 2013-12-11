#!/usr/bin/perl

=head1 NAME

I<rnaSeqDeNovoAssembly.pl>

=head1 DESCRIPTION

B<rnaSeqDeNovoAssembly.pl> Is the main de novo RNA assembly pipeline.

=head1 AUTHORS

B<David Morais> - I<dmorais@cs.bris.ac.uk>

B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

B<Joel Fillon> - I<joel.fillon@mcgill.ca>

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
use Cwd 'abs_path';
use Data::Dumper; # to delete
use Getopt::Std;
use LoadConfig;
use SampleSheet;
use SubmitToCluster;
use Trimmomatic;
use Trinity;

# Globals
#---------------------------------
my @A_steps = (
  {
    'name'   => 'trim',
    'loop'   => 'sample',
    'parent' => []
  },
  {
    'name'   => 'normalization',
    'loop'   => 'global',
    'parent' => ['trim']
  },
  {
    'name'   => 'trinity',
    'loop'   => 'global',
    'parent' => ['normalization']
  },
  {
    'name'   => 'blastSplitQuery',
    'loop'   => 'global',
    'parent' => ['trinity']
  },
  {
    'name'   => 'blast',
    'loop'   => 'global',
    'parent' => ['blastSplitQuery']
  },
  {
    'name'   => 'blastMergeResults',
    'loop'   => 'global',
    'parent' => ['blast']
  },
  {
    'name'   => 'rsemPrepareReference',
    'loop'   => 'global',
    'parent' => ['trinity']
  },
  {
    'name'   => 'rsem',
    'loop'   => 'sample',
    'parent' => ['rsemPrepareReference']
  },
  {
    'name'   => 'differentialGeneExpression',
    'loop'   => 'global',
    'parent' => ['rsem', 'blastMergeResults']
  }
);

# Create step hash indexed by step name for easy step retrieval
my %H_steps =  map {$_->{'name'} => $_} @A_steps;

my $configFile;
my $nanuqSampleSheet;
my $workDirectory;
my $designFile;

main();

# SUB
#-----------------
sub getUsage {
  my $usage = <<END;
Usage: perl $0 -h | -c CONFIG_FILE -s start_step_num -e end_step_num -n SAMPLE_SHEET [-d DESIGN_FILE] [-w WORK_DIR]
  -h  help and usage
  -c  .ini config file
  -s  start step, inclusive
  -e  end step, inclusive
  -n  nanuq sample sheet
  -d  design file
  -w  work directory (default current)

Steps:
END

  # List and number step names
  for (my $i = 1; $i <= @A_steps; $i++) {
    $usage .= $i . "- " . $A_steps[$i - 1]->{'name'} . "\n";
  }

  return $usage;
}

sub main {
  # Check options
  my %opts;
  getopts('hc:s:e:n:d:w:', \%opts);

  if (defined($opts{'h'}) ||
     !defined($opts{'c'}) ||
     !defined($opts{'s'}) ||
     !defined($opts{'e'}) ||
     !defined($opts{'n'})) {
    die (getUsage());
  }

  my $startStep = $opts{'s'};
  my $endStep = $opts{'e'};
  $configFile = $opts{'c'};
  $nanuqSampleSheet = $opts{'n'};
  $designFile = $opts{'d'};
  $workDirectory = $opts{'w'};

  # Get config and sample values
  my %cfg = LoadConfig->readConfigFile($configFile);
  my $rHoAoH_sampleInfo = SampleSheet::parseSampleSheetAsHash($nanuqSampleSheet);

  SubmitToCluster::initPipeline($workDirectory);

  for (my $i = $startStep; $i <= $endStep; $i++) {
    my $step = $A_steps[$i - 1];
    my $stepName = $step->{'name'};
    $step->{'jobIds'} = ();

    if ($step->{'loop'} eq 'sample') {
      foreach my $sample (keys %$rHoAoH_sampleInfo) {
        my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sample};
        &$stepName(\%cfg, $step, $workDirectory, $sample, $rAoH_sampleLanes);
      }
    } else {  # Global step
      &$stepName(\%cfg, $step, $workDirectory);
    }
  }
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

  # Retrieve the list of step parents
  my @A_stepParents = map {$H_steps{$_}} @{$step->{'parent'}};

  # Retrieve the list of lists of step parent job IDs if any
  my @AoA_stepParentJobIds = map {defined $_->{'jobIds'} ? $_->{'jobIds'} : []} @A_stepParents;

  # Flatten this list
  my @A_stepParentJobIds = map {@$_} @AoA_stepParentJobIds;

  # Concatenate all job IDs with cluster dependency separator
  $dependencies = join (getParam($rH_cfg, 'default', 'clusterDependencySep'), map {"\$" . $_} @A_stepParentJobIds);

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

sub trim {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;
  my $sample = shift;
  my $rAoH_sampleLanes = shift;

  for my $rH_laneInfo (@$rAoH_sampleLanes) {

    my $trimDirectory = "\$WORK_DIR/reads/$sample/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'};
    print 'mkdir -p ' . $trimDirectory . "\n";
    my $rO_job = Trimmomatic::trim($rH_cfg, $sample, $rH_laneInfo, $trimDirectory);

    submitJob($rH_cfg, $step, $sample, $rO_job);
  }
}

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

sub differentialGeneExpression {
  my $rH_cfg = shift;
  my $step = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();
  if (!$rO_job->isUp2Date()) {
    unless (defined $designFile) {die "Error: design file is not defined! (use -d option)\n" . getUsage()};
    unless (-f $designFile) {die "Error: design file $designFile does not exist!\n" . getUsage()};
    $designFile = abs_path($designFile);

    my $command = "\n";

    my $program = getParam($rH_cfg, 'blast', 'blastProgram');
    my $db = getParam($rH_cfg, 'blast', 'blastDb');
    my $blastDir = "\$WORK_DIR/blast";
    my $blastResult = "$blastDir/$program" . "_Trinity_$db.tsv";

    my $dgeDir = "\$WORK_DIR/DGE";
    my $isoformsMatrix = "$dgeDir/isoforms.counts.matrix";
    my $isoformsAnnotatedMatrix = "$dgeDir/isoforms.annotated.counts.matrix";
    my $genesMatrix = "$dgeDir/genes.counts.matrix";
    my $genesAnnotatedMatrix = "$dgeDir/genes.annotated.counts.matrix";

    $command .= moduleLoad($rH_cfg, [
      ['differentialGeneExpression', 'moduleVersion.trinity'],
      ['differentialGeneExpression', 'moduleVersion.cranR'],
      ['differentialGeneExpression', 'moduleVersion.tools']
    ]);

    $command .= "mkdir -p $dgeDir\n";

    # Create isoforms and genes matrices with counts of RNA-seq fragments per feature using Trinity RSEM utility
    $command .= "merge_RSEM_frag_counts_single_table.pl \$WORK_DIR/rsem/*/*.isoforms.results > $isoformsMatrix\n";
    $command .= "merge_RSEM_frag_counts_single_table.pl \$WORK_DIR/rsem/*/*.genes.results > $genesMatrix\n";

    # Merge isoforms and genes matrices with blast annotations if any
    $command .= "awk '!x[\\\$1]++' $blastResult | awk -F\\\"\\t\\\" 'FNR==NR {a[\\\$1]=\\\$2; next}{OFS=\\\"\\t\\\"; if (a[\\\$1]) {print \\\$1, a[\\\$1]} else {print \\\$1, \\\$1}}' - $isoformsMatrix | sed '1s/^\\t/Isoform\\tSymbol/' | paste - <(cut -f 2- $isoformsMatrix) | sed '1s/\\.isoforms\\.results//g' > $isoformsAnnotatedMatrix\n";
    $command .= "awk '!x[\\\$1]++' $blastResult | awk -F\\\"\\t\\\" 'FNR==NR {sub(/_seq.*/, \\\"\\\", \\\$1); a[\\\$1]=\\\$2; next}{OFS=\\\"\\t\\\"; if (a[\\\$1]) {print \\\$1, a[\\\$1]} else {print \\\$1, \\\$1}}' - $genesMatrix | sed '1s/^\\t/Gene\\tSymbol/' | paste - <(cut -f 2- $genesMatrix) | sed '1s/\\.genes\\.results//g' > $genesAnnotatedMatrix \n";

    $command .= "Rscript \\\$R_TOOLS/edger.R -d $designFile -c $isoformsAnnotatedMatrix -o $dgeDir \n";
    $command .= "Rscript \\\$R_TOOLS/edger.R -d $designFile -c $genesAnnotatedMatrix -o $dgeDir \\\n";

    $rO_job->addCommand($command);
  }
  submitJob($rH_cfg, $step, undef, $rO_job);
}
