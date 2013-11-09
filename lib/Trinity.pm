#!/usr/env/perl

=head1 NAME

I<Trinity>

=head1 SYNOPSIS

Trinity::sub(args)

B<Trinity::chrysalis>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $pair1, $pair2). In the case where there are multiple samples
the $pair1 and $pair2 will be  string with all the  samples (reads/sample1_left   reads/sample2_left).
 
B<Trinity::butterfly>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $fileButterflyComand)

B<Trinity::concatFastaCreateGtf>(%ref_hash_config, $sample_name, %ref_hash_laneInfo)


All subroutines return a ref_hash with the command line



=head1 DESCRIPTION

B<Trinity> is a library to use the
transcriptome assembly package, Trinity.

The lib has three subroutines: B<chrysalis, butterfly and concatFastaCreateGtf>.

- chrysalis and butterfly are steps of the package trinity.

- concatFastaCreateGtf is used to concatenate the fasta files and 
to generate de gtf file.

Each sub must be called independently with their set of variables.

=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug



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
use Data::Dumper;
use Config::Simple;
use LoadConfig;
use Cwd qw(abs_path);
use File::Basename;
use Job;

#-------------------
# SUB
#-------------------
sub normalize_by_kmer_coverage {
  my $rH_cfg = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();
#  $rO_job->testInputOutputs([$leftList, $rightList], ["$workDir/normalized_reads/both.fa"]);

  if (!$rO_job->isUp2Date()) {
    my $ram = "200G";
    my $CPU = "10";
    my $command;

    my $leftList = "\$WORK_DIR/reads/left_pair1.fastq.gz.list";
    my $rightList = "\$WORK_DIR/reads/right_pair2.fastq.gz.list";

    # Create sorted left/right lists of fastq.gz files
    $command .= "find \$WORK_DIR/reads/ -name *pair1*.fastq.gz | sort > $leftList; ";
    $command .= "find \$WORK_DIR/reads/ -name *pair2*.fastq.gz | sort > $rightList; ";

    # Load modules and run Trinity normalization
    $command .= 'module load ' . LoadConfig::getParam($rH_cfg, 'trinity', 'moduleVersion.trinity') . '; ';
    $command .= "normalize_by_kmer_coverage.pl \\
      --seqType fq \\
      --JM $ram \\
      --max_cov 30 \\
      --left_list $leftList \\
      --right_list $rightList \\
      --pairs_together \\
      --SS_lib_type RF \\
      --output \$WORK_DIR/normalization \\
      --JELLY_CPU $CPU \\
      --PARALLEL_STATS \\
      --KMER_SIZE 25 \\
      --max_pct_stdev 100";

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
    my $ram = "50G";
    my $CPU = 10;
    my $bflyCPU = int($CPU / 2);
    my $command;

    my $leftList = "\$WORK_DIR/reads/left_pair1.fastq.gz.list";
    my $rightList = "\$WORK_DIR/reads/right_pair2.fastq.gz.list";

    $command .= 'module load ' . LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.java') . ' ' .
      LoadConfig::getParam($rH_cfg, 'trinity', 'moduleVersion.trinity') . ' ' .
      LoadConfig::getParam($rH_cfg, 'bowtie', 'moduleVersion.bowtie') . ' ' .
      LoadConfig::getParam($rH_cfg, 'samtools', 'moduleVersion.samtools') . '; ';
    $command .= "Trinity.pl \\
      --seqType fq \\
      --JM $ram \\
      --left  $leftList.normalized_K25_C30_pctSD100.fq \\
      --right $rightList.normalized_K25_C30_pctSD100.fq \\
      --SS_lib_type RF \\
      --output \$WORK_DIR/trinity_out_dir \\
      --CPU $CPU \\
      --min_contig_length 200 \\
      --jaccard_clip \\
      --min_kmer_cov 2 \\
      --inchworm_cpu $CPU \\
      --bflyHeapSpaceMax 10G \\
      --bflyGCThreads 1 \\
      --bflyCPU $bflyCPU";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}


sub splitButterfly {
  my $rH_cfg  = shift;
  my $workDir = shift;
  my $butterflyCommandChunksNumber = shift;

  my $rO_job = new Job();
  my $butterflyCommands = "$workDir/trinity/chrysalis/butterfly_commands";

  $rO_job->testInputOutputs([$butterflyCommands], []);

  if (!$rO_job->isUp2Date()) {
    my $butterflyCommandsChunksDir = "$workDir/trinity/chrysalis/butterfly_commands_chunks";
    my $command = "mkdir -p $butterflyCommandsChunksDir; ";
    $command .= "shuf $butterflyCommands | split -d -a 3 -l \\`awk \'END{print int((NR - 1) / $butterflyCommandChunksNumber + 1)}\' $butterflyCommands\\` - $butterflyCommandsChunksDir/butterfly_commands.";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub rsemPrepareReference {
  my $rH_cfg = shift;
  my $workDirectory = shift;

  my $rO_job = new Job();

  if (!$rO_job->isUp2Date()) {
    my $command;

    $command .= 'module load ' . LoadConfig::getParam($rH_cfg, 'trinity', 'moduleVersion.trinity') . ' ' .
      LoadConfig::getParam($rH_cfg, 'bowtie', 'moduleVersion.bowtie') . ' ' .
      LoadConfig::getParam($rH_cfg, 'rsem', 'moduleVersion.rsem') . '; ';

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
    my $CPU = 4;
    my $command;

    my $left  = "\\`find \$WORK_DIR/reads -name $sample*pair1*.fastq.gz | sort | paste -s -d,\\`";
    my $right  = "\\`find \$WORK_DIR/reads -name $sample*pair2*.fastq.gz | sort | paste -s -d,\\`";

    $command .= 'module load ' . LoadConfig::getParam($rH_cfg, 'trinity', 'moduleVersion.trinity') . ' ' .
      LoadConfig::getParam($rH_cfg, 'bowtie', 'moduleVersion.bowtie') . ' ' .
      LoadConfig::getParam($rH_cfg, 'rsem', 'moduleVersion.rsem') . '; ';

    $command .= "run_RSEM_align_n_estimate.pl \\
      --transcripts \$WORK_DIR/trinity_out_dir/Trinity.fasta \\
      --left $left \\
      --right $right \\
      --seqType fq \\
      --SS_lib_type RF \\
      --prefix $sample \\
      --output_dir \$WORK_DIR/rsem/$sample \\
      --thread_count $CPU";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub mergeCounts {
  my $rH_cfg               = shift;
  my $rA_filePrefixToMerge = shift;
  my $outputIso            = shift;
  my $outputGene           = shift;

  my $rO_job = new Job();
  $rO_job->testInputOutputs($rA_filePrefixToMerge,[$outputIso, $outputGene]);
  
  if (!$rO_job->isUp2Date()) {
    my $command = undef;
    $command .= 'module load '.LoadConfig::getParam( $rH_cfg, 'mergeCounts', 'moduleVersion.trinity' ) . ' ;';
    $command .= ' \$TRINITY_HOME/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl ';
    for my $input (@{$rA_filePrefixToMerge}){
      $command .= ' '.$input.'.rsem.isoforms.results';
    }
    $command .= " | sed 's/\\.rsem\\.isoforms\\.results//g' > ".$outputIso;
    $command .= ' ; ';
    $command .= ' \$TRINITY_HOME/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl ';
    for my $input (@{$rA_filePrefixToMerge}){
      $command .= ' '.$input.'.rsem.genes.results';
    }
    $command .= " | sed 's/\\.rsem\\.genes\\.results//g' > ".$outputGene;


    $rO_job->addCommand($command);
  }

  return $rO_job;
}

1;
