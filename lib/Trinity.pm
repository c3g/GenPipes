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

  return "module load " . join(" ", map {getParam($rH_cfg, $_->[0], $_->[1])} @$rA_modules) . "\n";
}


sub normalize_by_kmer_coverage {
  my $rH_cfg = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();
#  $rO_job->testInputOutputs([$leftList, $rightList], ["$workDir/normalized_reads/both.fa"]);

  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    my $leftList = "\$WORK_DIR/reads/left_pair1.fastq.gz.list";
    my $rightList = "\$WORK_DIR/reads/right_pair2.fastq.gz.list";

    # Create sorted left/right lists of fastq.gz files
    $command .= "find \$WORK_DIR/reads/ -name *pair1*.fastq.gz | sort > $leftList\n";
    $command .= "find \$WORK_DIR/reads/ -name *pair2*.fastq.gz | sort > $rightList\n";

    # Load modules and run Trinity normalization
    $command .= moduleLoad($rH_cfg, [['trinity', 'moduleVersion.trinity']]);
    $command .= "normalize_by_kmer_coverage.pl \\
 --left_list $leftList \\
 --right_list $rightList \\
 --output \$WORK_DIR/normalization \\\n";
    $command .= " --JM " . getParam($rH_cfg, 'normalization', 'jellyfishMemory') . " \\\n";
    $command .= " --JELLY_CPU " . getParam($rH_cfg, 'normalization', 'jellyfishCPU') . " \\\n";
    $command .= " " . getParam($rH_cfg, 'normalization', 'normalizationOptions') . " \\\n";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub trinity {
  my $rH_cfg  = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();
#  $rO_job->testInputOutputs(["$workDir/normalized_reads/pairs.K25.stats.C30.pctSD100.accs"], []);

  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    my $leftList = "\$WORK_DIR/reads/left_pair1.fastq.gz.list";
    my $rightList = "\$WORK_DIR/reads/right_pair2.fastq.gz.list";

    $command .= moduleLoad($rH_cfg, [
      ['trinity', 'moduleVersion.java'],
      ['trinity', 'moduleVersion.trinity'],
      ['trinity', 'moduleVersion.bowtie'],
      ['trinity', 'moduleVersion.samtools']
    ]);

    $command .= "Trinity.pl \\
 --left  $leftList.normalized_*.fq \\
 --right $rightList.normalized_*.fq \\
 --output \$WORK_DIR/trinity_out_dir \\\n";
    $command .= " --JM " . getParam($rH_cfg, 'trinity', 'jellyfishMemory') . " \\\n";
    $command .= " --CPU " . getParam($rH_cfg, 'trinity', 'trinityCPU') . " \\\n";
    $command .= " --bflyCPU " . getParam($rH_cfg, 'trinity', 'bflyCPU') . " \\\n";
    $command .= " " . getParam($rH_cfg, 'trinity', 'trinityOptions') . " \\\n";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub trinityQC {
  my $rH_cfg  = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();
  if (!$rO_job->isUp2Date()) {
    my $command = "\n";

    $command .= moduleLoad($rH_cfg, [
      ['trinity', 'moduleVersion.cranR'],
      ['trinity', 'moduleVersion.trinity'],
      ['trinity', 'moduleVersion.bowtie'],
      ['trinity', 'moduleVersion.samtools']
    ]);

    $command .= "Rscript -e 'library(gqSeqUtils); dnaFastaStats(filename = \\\"\$WORK_DIR/trinity_out_dir/Trinity.fasta\\\", type = \\\"trinity\\\", output.prefix = \\\"\$WORK_DIR/trinity_out_dir/Trinity.stats\\\")' \\\n";
    #$command .= "alignReads.pl --seqType fa --left \$WORK_DIR/normalization/left.fa --right \$WORK_DIR/normalization/right.fa --SS_lib_type RF --retain_intermediate_files --aligner bowtie --target \$WORK_DIR/trinity_out_dir/Trinity.fasta -- -p 4 \\\n";

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
