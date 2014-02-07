#!/usr/env/perl

=head1 NAME

I<BVATools>

=head1 SYNOPSIS

BVATools->qc()

=head1 DESCRIPTION

B<BVATools> : Bam and Variant Analysis Tools


=head1 AUTHOR
B<Louis Letourneau> - I<louis.letourneau@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package BVATools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub qc {
  my $rH_cfg         = shift;
  my $read1          = shift;
  my $read2          = shift;
  my $type           = shift;
  my $regionName     = shift;
  my $outputDirectory= shift;
  
  my $ro_job = new Job();

  my $rA_inputs;
  if(defined($read2)) {
    $rA_inputs = [$read1, $read2];
  } else {
    $rA_inputs = [$read1];
  }
  $ro_job->testInputOutputs($rA_inputs, [$outputDirectory.'/mpsQC_'.$regionName.'_stats.xml']);
  
  if (!$ro_job->isUp2Date()) {
    my $nbThreads = LoadConfig::getParam($rH_cfg, 'generateQCGraphs','nbThreads');
    my $command;
    
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'generateQCGraphs','moduleVersion.java') . " " .  LoadConfig::getParam($rH_cfg, 'generateQCGraphs','moduleVersion.BVATools') . ' &&'; 
    $command .= ' java ' .LoadConfig::getParam($rH_cfg, 'generateQCGraphs', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'generateQCGraphs', 'maxRam').' -jar $BVATOOLS_JAR';
    $command .= ' readsqc --regionName \'' . $regionName . '\' --type ' . $type . ' --output \'' . $outputDirectory . '\' --read1 \'' . $read1 .'\'';
    if (defined($read2)) {
      $command .= ' --read2 \'' . $read2 . '\'';
    }
    if (defined($nbThreads)) {
      $command .= ' --threads ' . $nbThreads;
    }
    $ro_job->addCommand($command);
  }
  
  return $ro_job;
}



1;
