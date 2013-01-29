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
  my $outputIndexFile= $outputFolder. 'index.html';

  my $command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestFile) || !defined(-M $outputIndexFile) || $latestFile < -M $outputIndexFile) {
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'metrics','bwaModule') .' ;';
    $command .= ' module load ' .LoadConfig::getParam($rH_cfg, 'metrics','rnaseqModule') .' ;';
    $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'metrics', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'metrics', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'metrics', 'metricsRam').' -jar \${RNASEQC_JAR}';
    $command .= ' -n ' .LoadConfig::getParam($rH_cfg, 'metrics','topTranscript')
    $command .= ' -s ' .$inputFile;
    $command .= ' -t ' .LoadConfig::getParam($rH_cfg, 'metrics','referenceGtf');
    $command .= ' -r ' .LoadConfig::getParam($rH_cfg, 'metrics','referenceFasta');
    $command .= ' -o ' .$outputFolder ;
    $command .= ' -BWArRNA ' .LoadConfig::getParam($rH_cfg, 'metrics','ribosomalGtf');
  }
    
  return $command;
}

sub saturation {
	my $rH_cfg      = shift;
	my $countFile   = shift;
	my $gtfFile     = shift;
	my $rpkmDir = shift;
	my $saturationDir = shift;
	

	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'saturation' , 'cranRModule') .' ' . LoadConfig::getParam($rH_cfg, 'saturation' , 'toolsModule') . ' ;';
	$command .= ' Rscript $R_TOOLS/rpkmSaturation.R ' .$countFile .' ' .$gtfFile .' ' .$rpkmDir .' ' .$saturationDir;
	$command .= ' ' .LoadConfig::getParam($rH_cfg, 'saturation' , 'optionR');

	return $command;
}

sub fpkmCor {
	my $rH_cfg         = shift;
	my $paternFile     = shift;
	my $folderFile     = shift;
	my $outputBaseName = shift;
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'metrics' , 'cranRModule') .' ' . LoadConfig::getParam($rH_cfg, 'metrics' , 'toolsModule') . ' ;';
	$command .= ' Rscript $R_TOOLS/fpkmStats.R ' .$paternFile .' ' .$folderFile .' ' .$outputBaseName;
	
	return $command;
}

	

1;
