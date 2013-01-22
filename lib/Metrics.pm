#!/usr/env/perl

=head1 NAME

I<Metrics>

=head1 SYNOPSIS

Metrics-> rnaQc()

=head1 DESCRIPTION

B<Metrics> is a library to generate QC, stats and metrics

Input = file_name

Output = array


=head1 AUTHOR
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package Metrics;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------
sub rnaQc{
  my $rH_cfg        = shift;
  my $inputFile      = shift;
  my $outputFolder     = shift;


  my $latestFile = -M $inputFile;
  my $outputIndexFile= $outputFolder. 'index.html'

  my $command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestFile) || !defined(-M $outputIndexFile) || $latestFile < -M $outputIndexFile) {
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'metricsRNA','bwaModule') .' ;';
    $command .= ' module load ' .LoadConfig::getParam($rH_cfg, 'metricsRNA','rnaseqModule') .' ;';
    $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'metricsRNA', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'metricsRNA', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'metricsRNA', 'metricsRam').' -jar \${RNAQC_JAR}';
    $command .= ' -n ' .LoadConfig::getParam($rH_cfg, 'metricsRNA','topTranscript')
    $command .= ' -s ' .$inputFile;
    $command .= ' -t ' .LoadConfig::getParam($rH_cfg, 'metricsRNA','referenceGtf');
    $command .= ' -r ' .LoadConfig::getParam($rH_cfg, 'metricsRNA','referenceFasta');
    $command .= ' -o ' .$outputFolder ;
    $command .= ' -BWArRNA ' .LoadConfig::getParam($rH_cfg, 'metricsRNA','ribosomalGtf');
  }
    
  return $command;
}

 

1;
