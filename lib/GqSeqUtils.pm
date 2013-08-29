#!/usr/env/perl

=head1 NAME

I<GqSeqUtils>

=head1 SYNOPSIS

GqSeqUtils-> clientReport()

=head1 DESCRIPTION

B<GqSeqUtils> is a library to access/launch functions from the gqSeqUtils R package


=head1 AUTHOR
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package GqSeqUtils;

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
	my $pipeline = "";
	my $pipelineTMP = LoadConfig::getParam($rH_cfg, 'report','projectName');
	if (defined($pipelineTMP) && !($pipelineTMP eq "")) {
		$pipeline= 'pipeline=\"' .$pipelineTMP .'\",';
	}
	my $title = "";
	my $titleTMP = LoadConfig::getParam($rH_cfg, 'report','report.title');
	if (defined($titleTMP) && !($titleTMP eq "")) {
                $title= 'report.title=\"' .$titleTMP .'\",';
        }
	my $path = "";
	my $pathTMP = LoadConfig::getParam($rH_cfg, 'report','report.path');
        if (defined($pathTMP) && !($pathTMP eq "")) {
                $path= 'report.path=\"' .$pathTMP .'\",';
        }
	my $author = "";
	my $authorTMP = LoadConfig::getParam($rH_cfg, 'report','report.author');
	if (defined($authorTMP) && !($authorTMP eq "")) {
                $author= 'report.author=\"' .$authorTMP .'\",';
        }
        my $contact = "";
        my $contactTMP = LoadConfig::getParam($rH_cfg, 'report','report.contact');
        if (defined($contactTMP) && !($contactTMP eq "")) {
                $contact= 'report.contact=\"' .$contactTMP .'\",';
        }

	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'report','moduleVersion.cranR') .' &&';
	$command .= ' R --vanilla -e \'library(gqSeqUtils) ;';
	$command .= ' mugqicPipelineReport(';
	$command .= ' ' .$pipeline ;
	$command .= ' ' .$title ;
        $command .= ' ' .$path ;
        $command .= ' ' .$author ;
        $command .= ' ' .$contact ;
        $command .= ' ini.file.path=\"' .$iniFilePath .'\",' ;
        $command .= ' project.path=\"' .$projectPath .'\")\'' ;

	return $command 
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
