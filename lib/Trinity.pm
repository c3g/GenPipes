#!/usr/env/perl

=head1 NAME

I<Trinity>

=head1 SYNOPSIS

Trinity::sub(args)

=head1 DESCRIPTION

B<Trinity> is a library to use the
transcriptome assembly package Trinity.

=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk
Joel Fillon joel.fillon@mcgill.ca

=cut

package Trinity;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib $FindBin::Bin;

# Dependencies
#-----------------------
use LoadConfig;
use Job;

#-------------------
# SUB
#-------------------

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

  # Retrieve module values from the configuration file
  my @moduleValues = map {getParam($rH_cfg, $_->[0], $_->[1])} @$rA_modules;

  # Check by a system call if module is available
  for my $moduleValue (@moduleValues) {
    my $moduleShowOutput = `source /etc/profile.d/modules.sh; module show $moduleValue 2>&1`;
    $moduleShowOutput !~ /Error/i or die "Error in configuration file:\n$moduleShowOutput";
  }

  return "module load " . join(" ", map {getParam($rH_cfg, $_->[0], $_->[1])} @$rA_modules) . " && \\\n";
}


sub normalize_by_kmer_coverage {
  my $rH_cfg = shift;
  my $rA_leftReadFiles = shift;
  my $rA_rightReadFiles = shift;
  my $rA_singleReadFiles = shift;
  my $outputDirectory = shift;

  # Find out if reads are paired or single-end
  my $readType;
  if (defined($rA_leftReadFiles) and defined($rA_rightReadFiles) and not(defined($rA_singleReadFiles))) {
    $readType = "paired";
  } elsif (not(defined($rA_leftReadFiles)) and not(defined($rA_rightReadFiles)) and defined($rA_singleReadFiles)) {
    $readType = "single";
  } else {
    die "Error in normalize_by_kmer_coverage: mixed or undefined paired/single reads!\n";
  }

  my $maxCoverage = getParam($rH_cfg, 'normalization', 'maxCoverage');
  my $kmerSize = getParam($rH_cfg, 'normalization', 'kmerSize');
  my $maxPctStdev = getParam($rH_cfg, 'normalization', 'maxPctStdev');

  my $outputSuffix = ".normalized_K" . $kmerSize . "_C" . $maxCoverage . "_pctSD" . $maxPctStdev . ".fq";

  my $leftList;
  my $rightList;
  my $singleCat;
  my $rA_inputs;
  my $rA_outputs;
  my $readFileOptions;

  if ($readType eq "paired") {    # Paired reads
    $leftList = "$outputDirectory/left";
    $rightList = "$outputDirectory/right";

    $rA_inputs = [@$rA_leftReadFiles, @$rA_rightReadFiles];
    $rA_outputs = [$leftList . $outputSuffix, $rightList . $outputSuffix];
  } else {    # Single reads
    $singleCat = "$outputDirectory/single";
    $rA_inputs = [@$rA_singleReadFiles];
    $rA_outputs = [$singleCat . $outputSuffix];
  }

  my $rO_job = new Job();
  $rO_job->testInputOutputs($rA_inputs, $rA_outputs);

  if (!$rO_job->isUp2Date()) {
    my $command = "\n";
    $command .= "mkdir -p $outputDirectory && \\\n";

    # Create sorted left/right lists of fastq.gz files
    if ($readType eq "paired") {    # Paired reads
      # Check if same number of left and right reads
      @$rA_leftReadFiles == @$rA_rightReadFiles or die "Error in normalization: left and right files numbers differ!"; 

      $command .= "rm -f $leftList $rightList && \\\n";
      foreach my $leftReadFile (@$rA_leftReadFiles) {
        $command .= "echo $leftReadFile >> $leftList && \\\n";
      }
      foreach my $rightReadFile (@$rA_rightReadFiles) {
        $command .= "echo $rightReadFile >> $rightList && \\\n";
      }
      $readFileOptions = " --left_list $leftList --right_list $rightList ";
    } else {    # Single reads
      $command .= "rm -f $singleCat && \\\n";
      # Check if fastq are compressed or not
      my $catCmd;
      if ($$rA_singleReadFiles[0] =~ /\.gz$/) {$catCmd = "zcat"} else {$catCmd = "cat"};
      # Merge all single fastq in one file since trinityrnaseq_r20131110 does not support --single_list!
      $command .= "$catCmd " . join(" ", @$rA_singleReadFiles) . " > $singleCat && \\\n";
      $readFileOptions = " --single $singleCat ";
    }

    # Load modules and run Trinity normalization
    $command .= moduleLoad($rH_cfg, [['trinity', 'moduleVersion.trinity']]);
    $command .= "normalize_by_kmer_coverage.pl \\
$readFileOptions \\
 --output $outputDirectory \\\n";
    $command .= " --JM " . getParam($rH_cfg, 'normalization', 'jellyfishMemory') . " \\\n";
    $command .= " --JELLY_CPU " . getParam($rH_cfg, 'normalization', 'jellyfishCPU') . " \\\n";
    $command .= " --max_cov $maxCoverage \\\n";
    $command .= " --KMER_SIZE $kmerSize \\\n";
    $command .= " --max_pct_stdev $maxPctStdev \\\n";
    $command .= " " . getParam($rH_cfg, 'normalization', 'normalizationOptions') . " && \\\n";

    # Count normalized reads for stats
    $command .= "wc -l " . @$rA_outputs[0] . " | awk '{print \\\"# normalized $readType reads\\t\\\"\\\$1 / 4}' > $outputDirectory/normalization.stats.tsv \\\n";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub trinity {
  my $rH_cfg  = shift;
  my $rA_leftReadFiles = shift;
  my $rA_rightReadFiles = shift;
  my $rA_singleReadFiles = shift;
  my $outputDirectory = shift;

  defined($outputDirectory) or die "Error in trinity: outputDirectory is not defined!";

  # Find out if reads are paired or single-end
  my $readType;
  if (defined($rA_leftReadFiles) and defined($rA_rightReadFiles) and not(defined($rA_singleReadFiles))) {
    $readType = "paired";
  } elsif (not(defined($rA_leftReadFiles)) and not(defined($rA_rightReadFiles)) and defined($rA_singleReadFiles)) {
    $readType = "single";
  } else {
    die "Error in trinity: mixed or undefined paired/single reads!\n";
  }

  my $rA_inputs;
  my $readFileOptions;

  if ($readType eq "paired") {    # Paired reads
    $rA_inputs = [@$rA_leftReadFiles, @$rA_rightReadFiles];
    $readFileOptions = " --left " . join(" ", @$rA_leftReadFiles) . " --right " . join(" ", @$rA_rightReadFiles);
  } else {    # Single reads
    $rA_inputs = [@$rA_singleReadFiles];
    $readFileOptions = " --single " . join(" ", @$rA_singleReadFiles);
  }

  my $rO_job = new Job();
  $rO_job->testInputOutputs($rA_inputs, ["$outputDirectory/Trinity.fasta"]);

  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    $command .= moduleLoad($rH_cfg, [
      ['trinity', 'moduleVersion.java'],
      ['trinity', 'moduleVersion.trinity'],
      ['trinity', 'moduleVersion.bowtie'],
      ['trinity', 'moduleVersion.samtools'],
      ['trinity', 'moduleVersion.cranR']
    ]);

    $command .= "Trinity.pl \\
$readFileOptions \\
 --output $outputDirectory \\\n";
    $command .= " --JM " . getParam($rH_cfg, 'trinity', 'jellyfishMemory') . " \\\n";
    $command .= " --CPU " . getParam($rH_cfg, 'trinity', 'trinityCPU') . " \\\n";
    $command .= " --bflyCPU " . getParam($rH_cfg, 'trinity', 'bflyCPU') . " \\\n";
    $command .= " " . getParam($rH_cfg, 'trinity', 'trinityOptions') . " && \\\n";

    # Create Trinity FASTA ZIP file for future deliverables
    $command .= "gzip -c $outputDirectory/Trinity.fasta > $outputDirectory/Trinity.fasta.gz && \\\n";

    # Compute assembly stats
    $command .= "Rscript -e 'library(gqSeqUtils); dnaFastaStats(filename = \\\"$outputDirectory/Trinity.fasta\\\", type = \\\"trinity\\\", output.prefix = \\\"$outputDirectory/Trinity.stats\\\")' \\\n";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub rsemPrepareReference {
  my $rH_cfg = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();

  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    $command .= moduleLoad($rH_cfg, [
      ['trinity', 'moduleVersion.trinity'],
      ['trinity', 'moduleVersion.bowtie'],
      ['trinity', 'moduleVersion.rsem']
    ]);

    $command .= "run_RSEM_align_n_estimate.pl \\
      --transcripts \$WORK_DIR/trinity_out_dir/Trinity.fasta \\
      --just_prep_reference ";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub rsem {
  my $rH_cfg = shift;
  my $workDirectory = shift;
  my $sample = shift;

  my $rO_job = new Job();

  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    my $left  = "\\`find \$WORK_DIR/reads -name $sample*pair1*.fastq.gz | sort | paste -s -d,\\`";
    my $right  = "\\`find \$WORK_DIR/reads -name $sample*pair2*.fastq.gz | sort | paste -s -d,\\`";

    $command .= moduleLoad($rH_cfg, [
      ['trinity', 'moduleVersion.trinity'],
      ['trinity', 'moduleVersion.bowtie'],
      ['trinity', 'moduleVersion.rsem']
    ]);

    $command .= "run_RSEM_align_n_estimate.pl \\
      --transcripts \$WORK_DIR/trinity_out_dir/Trinity.fasta \\
      --left $left \\
      --right $right \\
      --seqType fq \\
      --SS_lib_type RF \\
      --prefix $sample \\
      --output_dir \$WORK_DIR/rsem/$sample \\
      --thread_count " . getParam($rH_cfg, 'rsem', 'rsemCPU') . " \\\n";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub edgeR {
  my $rH_cfg = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();

  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    $command .= moduleLoad($rH_cfg, [
      ['trinity', 'moduleVersion.trinity'],
      ['trinity', 'moduleVersion.cranR']
    ]);

    $command .= "mkdir -p \$WORK_DIR/edgeR; ";

    $command .= "merge_RSEM_frag_counts_single_table.pl \\
      \$WORK_DIR/rsem/*/*.genes.results \\
      > \$WORK_DIR/edgeR/genes.counts.matrix; ";

    $command .= "run_DE_analysis.pl \\
      --matrix \$WORK_DIR/edgeR/genes.counts.matrix \\
      --method edgeR \\
      --output \$WORK_DIR/edgeR ";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

1;
