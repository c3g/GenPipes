#!/usr/env/perl

=head1 NAME

I<BLAST>

=head1 SYNOPSIS

B<BLAST::align>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $FastaFile, $db)


B<BLAST::alignParallel>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $FastaFile, $db)

B<BLAST::bestHit>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $db)

I<B<The $db arg is optional, made to be used if more than one database was specified in the config file>>

B<BLAST::alignParallel> will depend on two scripts and their paths must be put on the config file

1- B<ParallelBlast>: a script written by David Morais 

2- B<fastasplit>: www.bcgsc.ca/downloads/parts/software/resources/src/exonerate-1.4.0/src/util/fastasplit.c

B<Both scripts need to be in the same dierctory.>

=head1 DESCRIPTION

B<BLAST> is a library to use the
alignment package, BLAST.

The lib implements two subroutines: one that performs BLAST with its on parallelization options;
and another that performs an external parallelization. 

In my benchmark the external parallelization was 25% faster but fell free to try it out yourself.


=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug.




=cut

package BLAST;

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
use File::Basename;
#-------------------
# SUB
#-------------------
our $rH_cfg;
our $sampleName;
our $rH_laneInfo;
our $fileFasta;

sub align {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $fileFasta   = shift;

    my $db = shift;

    # option used if more than one db was specified on the config file.
    # In this case $db should be passed as an argument
    #-----------------------------------------------------------------
    $rH_cfg->{'blast.db'} = defined($db) ? $db : $rH_cfg->{'blast.db'};

    my $outFile = $fileFasta;
    $outFile =~ s/\*_//;
    my $command = '';
    my %retVal;
    my $laneDirectory = "assembly/" . $sampleName . "/";
    my $input = $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'})  . '/' . $fileFasta;
    my $output = $laneDirectory . 'fasta_split/' . $outFile . '_BLASTOUT.txt';

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$input], [$output]);

    if (!$ro_job->isUp2Date()) {
      $command .= 'module add mugqic/blast/2.2.27+;';
      $command .= ' mkdir -p ' . $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'}) . ';';
      $command .= ' ' . $rH_cfg->{'blast.program'} . ' -num_threads ' . $rH_cfg->{'blast.nbThreads'};
      $command .= ' -query ' . $input;
      $command .= ' -db ' . $rH_cfg->{'blast.db'} . ' -out ' . $output;
      $command .= ' ' . $rH_cfg->{'blast.options'};

      $ro_job->addCommand($command);
    }
    return $ro_job;

}

sub alignParallel {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $fileFasta   = shift;

    my $db = shift;

    # option used if more than one db was specified on the config file.
    # In this case $db should be passed as an argument
    #------------------------------------------------------------------
    $rH_cfg->{'blast.db'} = defined($db) ? $db : $rH_cfg->{'blast.db'};
	
    my $outFile = $fileFasta;
    $outFile =~ s/\*_//;
    my $command = '';
    my %retVal;
    my $laneDirectory = "assembly/" . $sampleName . "/";
    my $input = $laneDirectory . 'fasta_split/' . $fileFasta;
    my $output = $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'})  . '/' . $outFile . '_BLASTOUT.txt';

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$input], [$output]);

    if (!$ro_job->isUp2Date()) {
      $command .= 'module add mugqic/blast/2.2.27+;';
      $command .= ' mkdir -p ' . $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'}) . ';';
      $command .= ' ' . $rH_cfg->{'blast.parallelBlast'} . ' --file ' . $input;
      $command .= ' --OUT ' . $output;
      $command .= ' -n ' . $rH_cfg->{'blast.nbThreads'} . ' --BLAST ' . '\'' . $rH_cfg->{'blast.program'};
      $command .= ' -db ' . $rH_cfg->{'blast.db'} . ' ' . $rH_cfg->{'blast.options'} . '\'';

      $ro_job->addCommand($command);
    }
    return $ro_job;

}

sub bestHit {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;

    my $db = shift;

    # option used if more than one db was specified on the config file.
    # In this case $db should be passed as an argument
    #------------------------------------------------------------------
    $rH_cfg->{'blast.db'} = defined($db) ? $db : $rH_cfg->{'blast.db'};
    my $laneDirectory = "assembly/" . $sampleName . "/";
    my $input = $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'})  . '/';
    my $output = $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'})  . '/blast_BestHit.txt';

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$input], [$output]);

    if (!$ro_job->isUp2Date()) {
      my $command .= 'cat ' . $input . '*.txt ';
      $command .= ' >' . $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'})  . '/blastRes.txt ;';
      $command .= ' sh ' . $rH_cfg->{'blast.blastHq'} . ' ' . $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'})  . '/blastRes.txt ';
      $command .= ' >' . $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'})  . '/blastRes_HQ.txt ;';
      $command .= ' ' . $rH_cfg->{'blast.BestHit'} . ' ' . $laneDirectory . 'fasta_split/' . basename($rH_cfg->{'blast.db'})  . '/blastRes_HQ.txt ';
      $command .= ' >' . $output;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub dcmegablast{ #JT: Initially for for PacBio pipeline
	my $rH_cfg      = shift;
	my $gunzipCmd   = shift;
	my $infileFasta = shift; 
	my $outfmt      = shift;
	my $outfile     = shift;

	my $ro_job = new Job();
    $ro_job->testInputOutputs([$infileFasta], [$outfile]);
    
	if (!$ro_job->isUp2Date()) {
		my $cmd;
		
		$cmd .= 'module load '.LoadConfig::getParam($rH_cfg, 'blast', 'moduleVersion.blast').' ; ';
		$cmd .= $gunzipCmd . ' &&';
		$cmd .= ' blastn';
		$cmd .= ' -task dc-megablast';
		$cmd .= ' -query ' . $infileFasta;
		$cmd .= ' -outfmt ' . $outfmt;
		$cmd .= ' -out ' .$outfile;
		$cmd .= ' -num_threads ' .LoadConfig::getParam($rH_cfg, 'blast', 'num_threads');
		$cmd .= ' -db ' .LoadConfig::getParam($rH_cfg, 'blast', 'blastdb');


      $ro_job->addCommand($cmd);
    }
    return $ro_job;
}

sub blastdbcmd{ # JT: Initially for PacBio pipeline
	my $rH_cfg      = shift;
	my $entryCmd    = shift;
	my $outfile     = shift;

	my $ro_job = new Job();
    $ro_job->testInputOutputs(undef, [$outfile]);
    
	if (!$ro_job->isUp2Date()) {
		my $cmd;
		
		$cmd .= 'module load '.LoadConfig::getParam($rH_cfg, 'blast', 'moduleVersion.blast').' ; ';
		$cmd .= ' blastdbcmd';
		$cmd .= ' -db ' . LoadConfig::getParam($rH_cfg, 'blast', 'blastdb');
		$cmd .= ' -entry ' . $entryCmd;
		$cmd .= ' -outfmt %f';
		$cmd .= ' > ' . $outfile;

      $ro_job->addCommand($cmd);
    }
    return $ro_job;
}

1;

