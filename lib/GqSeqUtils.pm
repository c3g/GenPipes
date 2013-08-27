#!/usr/env/perl

=head1 NAME

I<GqSeqUtils>

=head1 SYNOPSIS

GqSeqUtils-> rnaReport()

=head1 DESCRIPTION

B<GqSeqUtils> is a library to access/launch functions from the gqSeqUtils R package


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
sub clientReport{
	my $rH_cfg        = shift;
	my $iniFilePath   = shift;
	my $projectPath   = shift;
	
# 	pipeline=        ini.file.path=   report.title=    report.contact=  
# project.path=    report.path=     report.author=   report.css=
	my $pipeline = " ";
	my $pipelineTMP = LoadConfig::getParam($rH_cfg, 'report','projectName');
	if (!defined($pipelineTMP) || $pipelineTMP eq "") {
		$pipeline= 'pipeline=' .$pipelineTMP;
	}
	my $title = LoadConfig::getParam($rH_cfg, 'report','report.title');
	my $path = LoadConfig::getParam($rH_cfg, 'report','report.path');
	my $author = LoadConfig::getParam($rH_cfg, 'report','report.author');
	
	

	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'report','moduleVersion.cranR') .' &&';
	$command .= ' R --vanilla -e \"library(gqSeqUtils) ;"';
	$command .= '' 
}


sub exploratoryRnaAnalysis{
        my $rH_cfg        = shift;
        my $readSetSheet  = shift;
	my $workDirectory = shift;
        my $configFile    = shift;

	my $rscript = 'suppressPackageStartupMessages(library(gqSeqUtils));';
        $rscript .= ' initIllmSeqProject(nanuq.file= \"' . $readSetSheet . '\",overwrite.sheets=FALSE,project.path= \"' . $workDirectory . '\" );';
        $rscript .= ' exploratoryRNAseq(project.path= \"' . $workDirectory . '\",ini.file.path = \"' . $configFile . '\" );';
        $rscript .= ' print(\"done.\")';
        my $command = 'module load ' .LoadConfig::getParam($rH_cfg, 'downstreamAnalyses','moduleVersion.cranR') .' &&';
	$command .= ' Rscript -e ' . '\'' . $rscript .'\'';

	return $command;
}

1;
