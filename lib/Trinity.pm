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

#-------------------
# SUB
#-------------------
sub chrysalis {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;
    my $pair1       = shift;    # For single command the left will receive the file.
    my $pair2       = shift;

    my $ro_job;
    if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
        $ro_job = _chrysalisSingleCommand($rH_cfg, $sampleName, $rH_laneInfo, $pair1);
    }
    elsif ( $rH_laneInfo->{'runType'} eq "PAIRED_END" ) {
        $ro_job = _chrysalisPairCommand($rH_cfg, $sampleName, $rH_laneInfo, $pair1, $pair2);
    }
    else {
        die "Unknown runType: " . $rH_laneInfo->{' runType '} . "\n";
    }
    return $ro_job;
}

sub butterfly {
    my $rH_cfg                = shift;
    my $sampleName            = shift;
    my $rH_laneInfo           = shift;
    my $fileButterflyCommand  = shift;

    my $laneDirectory = "assembly/" . $sampleName . "/chrysalis/";

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$laneDirectory . 'butterfly_split/' . $fileButterflyCommand], undef);

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
      $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
      $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
      $command .= ' ' . $rH_cfg->{'butterfly.parallel'} . ' -f ' . $laneDirectory . 'butterfly_split/' . $fileButterflyCommand;
      $command .= ' -n ' . $rH_cfg->{'butterfly.nbThreads'} . ' ';

      $ro_job->addCommand($command);
    }
    return $ro_job;

}

sub concatFastaCreateGtf {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;

    my $laneDirectory = "assembly/" . $sampleName . "/";

    my $ro_job = new Job();
    $ro_job->testInputOutputs(undef, [$laneDirectory . 'Trinity.fasta', $laneDirectory . 'Trinity.2.fasta']);

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
      $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
      $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
      $command .= ' find ' . $laneDirectory . 'chrysalis';
      $command .= ' -name "*allProbPaths.fasta" -exec cat {} + >' . $laneDirectory . 'Trinity.fasta &&';
      $command .= ' sh ' . $rH_cfg->{'trinity.createGtf'} . ' ' . $laneDirectory . 'Trinity.fasta';
      $command .= ' ' . $laneDirectory . $sampleName . '.gtf &&';
      $command .= ' awk \'{print \$1} \' ' . $laneDirectory . 'Trinity.fasta ';
      $command .= ' >' . $laneDirectory . 'Trinity.2.fasta';

      $ro_job->addCommand($command);
    }
    return $ro_job;

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

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$pair1, $pair2],[$outputFile]);
    
    if (!$ro_job->isUp2Date()) {
      $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
      $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
      $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
      $command .= ' Trinity.pl --seqType fq --JM 100G';
      $command .= ' --left' . ' \" ' . $pair1 . ' \" ' . '--right' . ' \" ' . $pair2 . ' \" ';
      $command .= ' --CPU ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'nbThreads');
      $command .= ' --output ' . $laneDirectory;
      $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

      $ro_job->addCommand($command);
    }
    return $ro_job;
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
    my $ro_job = new Job();
    $ro_job->testInputOutputs([$pair1],[$outputFile]);
    
    if (!$ro_job->isUp2Date()) {
      $command .= 'module add ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.java' );
      $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.bowtie' );
      $command .= ' ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'moduleVersion.trinity' ) . ' ;';
      $command .= ' Trinity.pl --seqType fq --JM 100G';
      $command .= ' ' . $pair1;
      $command .= ' --CPU ' . LoadConfig::getParam( $rH_cfg, 'trinity', 'nbThreads');
      $command .= ' --output ' . $laneDirectory;
      $command .= ' --min_kmer_cov 31 --max_reads_per_loop 200000000 --no_run_butterfly ';

      $ro_job->addCommand($command);
    }

    return $ro_job;

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

    my $ro_job = new Job();
    $ro_job->testInputOutputs([\@inputs],[$assembly]);
    
    if (!$ro_job->isUp2Date()) {
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

      $ro_job->addCommand($command);
    }

    return $ro_job;
}

sub mergeCounts {
    my $rH_cfg               = shift;
    my $rA_filePrefixToMerge = shift;
    my $outputIso            = shift;
    my $outputGene           = shift;

    my $ro_job = new Job();
    $ro_job->testInputOutputs($rA_filePrefixToMerge,[$outputIso, $outputGene]);
    
    if (!$ro_job->isUp2Date()) {
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


      $ro_job->addCommand($command);
    }

    return $ro_job;
}

1;
