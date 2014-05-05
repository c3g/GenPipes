#!/usr/bin/env perl

=head1 NAME

I<Cufflinks>

=head1 SYNOPSIS

Cleaning-> rna()

=head1 DESCRIPTION

B<Cleaning> is a library to do: pipeline temporary file cleaning

Input = file_name

Output = array


=head1 AUTHOR
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package Cleaning;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin";

# Dependencies
#-----------------------
use LoadConfig;
use File::Path;

# SUB
#-----------------------
sub rna {
  ### clean current directory for RNAseq tempory files
  &rawReads;
  &reads;
  &tophat;
  &fpkm;
  &cuffdiff;
  &raw_counts;
  &deliverable;
  &exploratory;
  &metrics;
  &tracks;
}

sub rawReads {
  ### cleaning raw_reads if fastq are only symlink
  my $rrClean = 1 ;
  foreach my $name ( <raw_reads/*/*/*.gz> ) {
    $rrClean = 0  if( ! -l $name );
  }
  if ( $rrClean ) {
    rmtree('raw_reads');
    print "the raw_reads folder is now removed\n";
  }
  else {
    print "the raw_reads folder contains real fastq !! Please clean it manually and carefully\n";
  }
  print "----------------------------------\n";
}

sub reads {
  print "Cleaning the reads folder:\n";
  foreach my $name ( <reads/*/*/*.gz> ) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <reads/*/output_jobs>) {
    print "$name \n";
    rmtree($name);
  }
  print "\n----------------------------------\n";
  print "the reads folder is now cleaned\n";
  print "----------------------------------\n";
}

sub tophat {
  print "Cleaning the alignment folder:\n";
  foreach my $name ( <alignment/*/*/*> ) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <alignment/*/*karyotypic.ba[im]>) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <alignment/*/*queryNameSorted.ba[im]>) {
    print "$name \n";
    rmtree($name);
  }
  print "\n----------------------------------\n";
  print "the alignment folder is now cleaned\n";
  print "----------------------------------\n";
}

sub fpkm {
  print "Cleaning the fpkm folder:\n";
  foreach my $name ( <fpkm/*/output_jobs> ) {
    print "$name \n";
    rmtree($name);
  }
  print "fpkm/output_jobs\n";
  rmtree("fpkm/output_jobs");
  foreach my $name ( <fpkm/*/*/*[gfa]>) {
    print "$name \n";
    rmtree($name);
  }
  print "\n----------------------------------\n";
  print "the fpkm folder is now cleaned\n";
  print "----------------------------------\n";
}

sub cuffdiff {
  print "Cleaning the cuffdiff folder:\n";
  foreach my $name ( <cuffdiff/*/output_jobs> ) {
    print "$name \n";
    rmtree($name);
  }
  print "cuffdiff/output_jobs\n";
  rmtree("cuffdiff/output_jobs");
  foreach my $name ( <cuffdiff/*/*/*tracking> ) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <cuffdiff/*/*/*info> ) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <cuffdiff/*/*/*diff> ) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <cuffdiff/*/*/*gtf> ) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <cuffdiff/*/*/*list> ) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <cuffdiff/*/*/logs> ) {
    print "$name \n";
    rmtree($name);
  }
  print "\n----------------------------------\n";
  print "the cuffdiff folder is now cleaned\n";
  print "----------------------------------\n";
}

sub raw_counts {
  rmtree('raw_counts');
  print "the raw_counts folder is now removed\n";
}

sub dge {
  print "Cleaning the DGE folder:\n";
  foreach my $name ( <DGE/*/output_jobs> ) {
    print "$name \n";
    rmtree($name);
  }
  print "DGE/output_jobs\n";
  rmtree("DGE/output_jobs");
  foreach my $name ( <DGE/*/deseq_results.csv>) {
    print "$name \n";
    rmtree($name);
  }
  foreach my $name ( <DGE/*/edger_results.csv>) {
    print "$name \n";
    rmtree($name);
  }
  print "\n----------------------------------\n";
  print "the DGE folder is now cleaned\n";
  print "----------------------------------\n";
}

sub deliverable {
  rmtree('deliverable');
  print "the deliverable folder is now removed\n";
}

sub exploratory {
  print "exploratory/output_jobs\n";
  rmtree("exploratory/output_jobs");
  print "\n----------------------------------\n";
  print "the exploratory folder is now cleaned\n";
  print "----------------------------------\n";
}

sub metrics {
  print "Cleaning the metrics folder:\n";
  foreach my $name ( grep {!/zip\|stats/} <raw_reads/*> ) {
    print "$name \n";
    rmtree($name);
  }
   print "\n----------------------------------\n";
  print "the exploratory folder is now cleaned\n";
  print "----------------------------------\n";
}
  
sub tracks {
  rmtree('tracks');
  print "the tracks folder is now removed\n";
}




sub cuffcompare {
  my $rH_cfg           = shift;
  my $rA_mergeList     = shift;
  my $outputPrefix     = shift;
  my $mergeGtfFilePath = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs($rA_mergeList, [$outputPrefix . '.combined.gtf', $outputPrefix . '.TranscriptList.tsv']);

  if (!$ro_job->isUp2Date()) {
    my $mergeListString = join(' ', @{$rA_mergeList});
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['cuffcompare', 'moduleVersion.cufflinks'], ['cuffcompare', 'moduleVersion.tools']]) . ' &&';
    $command .= ' cuffcompare -o ' . $outputPrefix;
    $command .= ' -r ' . LoadConfig::getParam($rH_cfg, 'cuffcompare', 'referenceGtf', 1, 'filepath');
    $command .= ' -R ' . LoadConfig::getParam($rH_cfg, 'fpkm', 'referenceFasta', 1, 'filepath');
    $command .= ' -T ' . $mergeListString . ' &&';
    $command .= ' formatDenovoCombinedGTF.py' ;
    $command .= ' -c ' . $outputPrefix . '.combined.gtf';
    $command .= ' -t ' . $outputPrefix . '.tracking';
    $command .= ' -s ' . $mergeGtfFilePath;
    $command .= ' -o ' . $outputPrefix . '.TranscriptList.tsv';

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

#sub cuffmerge {
#  my $rH_cfg        = shift;
#  my $mergeListFile = shift;
#  my $outputDir     = shift;
#  
#  my $ro_job = new Job();
#  $ro_job->testInputOutputs([], []);
#
#  if (!$ro_job->isUp2Date()) {
#    my $command;
#    $command .= LoadConfig::moduleLoad($rH_cfg, [['cuffmerge', 'moduleVersion.cufflinks']]) . ';';
#    $command .= ' cuffmerge -p ' . LoadConfig::getParam($rH_cfg, 'cuffmerge', 'numThreads', 1, 'int');
#    $command .= ' -o ' . $outputDir;
#    $command .= ' -g ' . LoadConfig::getParam($rH_cfg, 'cuffmerge', 'referenceGtf', 1, 'filepath');
#    $command .= ' -s ' . LoadConfig::getParam($rH_cfg, 'fpkm', 'referenceFasta', 1, 'filepath');
#    $command .= ' ' . $mergeListFile;
#  }
#  return $command;
#}


sub mergeGtfFormat {
  my $rH_cfg     = shift;
  my $inputFile  = shift;
  my $outputFile = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputFile], [$outputFile]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['cuffmerge', 'moduleVersion.tools']]) . ' &&';
    $command .= ' perl \$PERL_TOOLS/formatGtfCufflinks.pl' . ' '. $inputFile . ' ' . $outputFile ;

    $ro_job->addCommand($command);
  }

  return $ro_job;
}

sub mergeCuffdiffRes {
  my $rH_cfg     = shift;
  my $designFile = shift;
  my $outputDir  = shift;
  my $inputDir   = shift;

  my $ro_job = new Job();
  # Can't test directories!
  #$ro_job->testInputOutputs([$inputDir], [$outputDir]);

  if (!$ro_job->isUp2Date()) {
    ### TO DO: re-write mergecuffdiff_known.R and mergecuffdiff_denovo.R to be more portable
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['cuffmerge', 'moduleVersion.tools'], ['cuffdiff', 'moduleVersion.cranR']]) . ' &&';
    $command .= ' Rscript \$R_TOOLS/mergecuffdiff_known.R ' . $outputDir . ' ' . $inputDir . ' ' . $designFile ;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub filterResults {
  my $rH_cfg    = shift;
  my $outputDir = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs(undef, undef);

  if (!$ro_job->isUp2Date()) {
    ### TO DO : make it more portable when the mergeCuffdiffRes R script will be re-write
    my $command;
    $command .= 'for i in \`ls ' . $outputDir . '/*/isoform_exp.diff.with.fpkm.csv\` ;';
    $command .= ' do head -1 \$i > ' . $outputDir . '/tmp ;';
    $command .= ' sed 1d \$i | grep -v \"NOTEST\" | grep -v \"FAIL\" | sort -k 12 -g >> ' . $outputDir . '/tmp ;';
    $command .= ' mv ' . $outputDir . '/tmp \$i ;';
    $command .= ' done';

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

1;
