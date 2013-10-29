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
  my $rH_cfg    = shift;
  my $workDir   = shift;
  my $leftList  = shift;    # For single command the left will receive the file.
  my $rightList = shift;

  my $rO_job = new Job();
  $rO_job->testInputOutputs([$leftList, $rightList], ["$workDir/normalized_reads/both.fa"]);

  if (!$rO_job->isUp2Date()) {
    my $ram = "200G";
    my $ncores = "10";
    my $command;
    $command .= 'module load ' . LoadConfig::getParam($rH_cfg, 'trinity', 'moduleVersion.trinity') . '; ';
    $command .= "normalize_by_kmer_coverage.pl \\
      --seqType fq \\
      --JM $ram \\
      --max_cov 30 \\
      --left_list $leftList \\
      --right_list $rightList \\
      --pairs_together \\
      --SS_lib_type RF \\
      --output normalized_reads \\
      --JELLY_CPU $ncores \\
      --PARALLEL_STATS \\
      --KMER_SIZE 25 \\
      --max_pct_stdev 100";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub trinity {
  my $rH_cfg  = shift;
  my $workDir = shift;
  my $leftList  = shift;    # For single command the left will receive the file.
  my $rightList = shift;

  my $rO_job = new Job();
  $rO_job->testInputOutputs(["$workDir/normalized_reads/pairs.K25.stats.C30.pctSD100.accs"], []);

  if (!$rO_job->isUp2Date()) {
    my $ram = "500G";
    my $ncores = "20";
    my $command;
    $command .= 'module load ' . LoadConfig::getParam($rH_cfg, 'trinity', 'moduleVersion.trinity') . ' ' .
      LoadConfig::getParam($rH_cfg, 'bowtie', 'moduleVersion.bowtie') . ' ' .
      LoadConfig::getParam($rH_cfg, 'samtools', 'moduleVersion.samtools') . '; ';
    $command .= "Trinity.pl \\
      --seqType fq \\
      --JM $ram \\
      --left   $leftList.normalized_K25_C30_pctSD100.fq \\
      --right  $rightList.normalized_K25_C30_pctSD100.fq \\
      --SS_lib_type RF \\
      --output $workDir/assembly \\
      --CPU $ncores \\
      --min_contig_length 200 \\
      --jaccard_clip \\
      --min_kmer_cov 2 \\
      --inchworm_cpu $ncores \\
      --bflyCPU $ncores \\
      --no_run_butterfly";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub chrysalis {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rH_laneInfo = shift;
  my $pair1       = shift;    # For single command the left will receive the file.
  my $pair2       = shift;

  my $rO_job;
  if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
      $rO_job = _chrysalisSingleCommand($rH_cfg, $sampleName, $rH_laneInfo, $pair1);
  }
  elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
      $rO_job = _chrysalisPairCommand($rH_cfg, $sampleName, $rH_laneInfo, $pair1, $pair2);
  }
  else {
      die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
  }
  return $rO_job;
}

sub splitButterfly {
  my $rH_cfg  = shift;
  my $workDir = shift;
  my $butterflyCommandChunksNumber = shift;

  my $rO_job = new Job();
  my $butterflyCommands = "$workDir/assembly/chrysalis/butterfly_commands";

  $rO_job->testInputOutputs([$butterflyCommands], []);

  if (!$rO_job->isUp2Date()) {
    my $butterflyCommandsChunksDir = "$workDir/assembly/chrysalis/butterfly_commands_chunks";
    my $command = "mkdir -p $butterflyCommandsChunksDir; ";
    $command .= "split -d -a 3 -l \\`awk \'END{print int((NR - 1) / $butterflyCommandChunksNumber + 1)}\' $butterflyCommands\\` $butterflyCommands $butterflyCommandsChunksDir/butterfly_commands.";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub butterfly {
  my $rH_cfg                 = shift;
  my $butterflyCommandsChunk = shift;

  my $rO_job = new Job();
  #$rO_job->testInputOutputs([$laneDirectory . 'butterfly_commands_chunks/' . $fileButterflyCommand], undef);

  if (!$rO_job->isUp2Date()) {
    my $command;
    $command .= 'module add ' . LoadConfig::getParam($rH_cfg, 'trinity', 'moduleVersion.java');
    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'trinity', 'moduleVersion.trinity') . ' ;';
    $command .= " if [ -f $butterflyCommandsChunk ]; then sh $butterflyCommandsChunk ; fi";

    $rO_job->addCommand($command);
  }
  return $rO_job;

}

sub concatFastaCreateGtf {
  my $rH_cfg      = shift;
  my $workDir     = shift;
#  my $sampleName  = shift;
#  my $rH_laneInfo = shift;

#  my $laneDirectory = "assembly/" . $sampleName . "/";

  my $rO_job = new Job();
  $rO_job->testInputOutputs(undef, [$workDir . '/assembly/Trinity.fasta', $workDir . '/assembly/Trinity.2.fasta']);

  if (!$rO_job->isUp2Date()) {
    my $assemblyDir = "$workDir/assembly";
    my $command;
    # Merge Butterfly output files into final Trinity fasta assembly
    $command .= "find $assemblyDir/chrysalis -name *allProbPaths.fasta -exec cat {} + > $assemblyDir/Trinity.fasta && ";

    # Convert fasta to gtf
    $command .= "grep '^>' $assemblyDir/Trinity.fasta | perl -pe 's/^>((\\S+)_seq\\d+)\\s+len=(\\d+).*/\1\tprotein_coding\texon\t1\t\3\t.\t+\t.\tgene_id \"\2\"; transcript_id \"\1\";/' > $assemblyDir/Trinity.gtf && ";

    # Create fasta with no description
    $command .= "awk '{print \$1}' $assemblyDir/Trinity.fasta > $assemblyDir/Trinity.no_desc.fasta";

    $rO_job->addCommand($command);
  }
  return $rO_job;

}

sub _chrysalisPairCommand {
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rH_laneInfo = shift;
  my $pair1 = shift;
  my $pair2 = shift;

  my $command = '';
  my %retVal;

  my $laneDirectory = 'assembly/' . $sampleName . '/';
  my $outputFile = $laneDirectory .'/chrysalis/butterfly_commands.adj';

  my $rO_job = new Job();
  $rO_job->testInputOutputs([$pair1, $pair2],[$outputFile]);
  
  if (!$rO_job->isUp2Date()) {
    $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
    $command .= ' Trinity.pl --seqType fq --JM 100G';
    $command .= ' --left' . ' \" ' . $pair1 . ' \" ' . '--right' . ' \" ' . $pair2 . ' \" ';
    $command .= ' --CPU ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'nbThreads');
    $command .= ' --output ' . $laneDirectory;
    $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub _chrysalisSingleCommand {
  my $rH_cfg = shift;
  my $sampleName = shift;
  my $rH_laneInfo = shift;
  my $pair1 = shift;

  my $command = '';
  my %retVal;
  
  my $laneDirectory = 'assembly/' . $sampleName . '/';
  my $outputFile = $laneDirectory .'/chrysalis/butterfly_commands.adj';
  my $rO_job = new Job();
  $rO_job->testInputOutputs([$pair1],[$outputFile]);
  
  if (!$rO_job->isUp2Date()) {
    $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
    $command .= ' Trinity.pl --seqType fq --JM 100G';
    $command .= ' ' . $pair1;
    $command .= ' --CPU ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'nbThreads');
    $command .= ' --output ' . $laneDirectory;
    $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

    $rO_job->addCommand($command);
  }

  return $rO_job;

}

sub abundance {
  my $rH_cfg       = shift;
  my $assembly     = shift;
  my $outputPrefix = shift;
  my $pair1        = shift;    # For single command the left will receive the file.
  my $pair2        = shift;

  my $unzippedPair1;
  my $unzippedPair2;
  my @inputs;

  $outputPrefix = abs_path($outputPrefix);
  $pair1 = abs_path($pair1);
  $pair1 =~ /(.+)\.gz/;
  $unzippedPair1 = $1;
  push(@inputs, $pair1);
  
  if(defined($pair2)) {
    $pair2 = abs_path($pair2);
    $pair2 =~ /(.+)\.gz/;
    $unzippedPair2 = $1;
    push(@inputs, $pair2);
  }

  my $outputFile = $outputPrefix.'.transcript.bam';

  my $rO_job = new Job();
  $rO_job->testInputOutputs([\@inputs],[$assembly]);
  
  if (!$rO_job->isUp2Date()) {
    my $command;
    $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'abundance', 'moduleVersion.java' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'abundance', 'moduleVersion.bowtie' );
    $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'abundance', 'moduleVersion.trinity' ) . ' ;';
    $command .= ' cd '.dirname($assembly).' ;';
    # We tried with mkfifo but bowtie keeps giving Broken Pipes
    $command .= ' gunzip -c '.$pair1.' > '.$unzippedPair1.' ;';
    if(defined($pair2)) {
      $command .= ' gunzip -c '.$pair2.' > '.$unzippedPair2.' ;';
    }
    $command .= ' \$TRINITY_HOME/util/RSEM_util/run_RSEM_align_n_estimate.pl ';
    $command .= ' --transcripts ' . basename($assembly);
    $command .= ' --seqType fq';
    if(!defined($pair2)) {
      $command .= ' --single '.$unzippedPair1;
    }
    else {
      $command .= ' --left '.$unzippedPair1;
      $command .= ' --right '.$unzippedPair2;
    }
    $command .= ' --thread_count '. LoadConfig::getParam( $rH_cfg, 'abundance', 'nbThreads' );
    $command .= ' --prefix ' . $outputPrefix;
    $command .= ' -- --bowtie-chunkmbs ' . LoadConfig::getParam( $rH_cfg, 'abundance', 'chunkmbs' ).';';
    $command .= ' rm '.$unzippedPair1.' ;';
    if(defined($pair2)) {
      $command .= ' rm '.$unzippedPair2.' ;';
    }

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
