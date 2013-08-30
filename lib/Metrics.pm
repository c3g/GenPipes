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
use PipelineUtils;
use LoadConfig;
use SAMtools;

# SUB
#-----------------------
sub rnaQc{
  my $rH_cfg        = shift;
  my $inputFile      = shift;
  my $outputFolder     = shift;


  my $outputIndexFile= $outputFolder. 'index.html';

  my $up2date = PipelineUtils::testInputOutputs([$inputFile], [$outputFolder .'.zip', $outputIndexFile]);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'rnaQc','moduleVersion.bwa') ;
    $command .= ' ' .LoadConfig::getParam($rH_cfg, 'rnaQc','moduleVersion.rnaseqc') .' &&';
    $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'rnaQc', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'rnaQc', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'rnaQc', 'metricsRam').' -jar \${RNASEQC_JAR}';
    $command .= ' -n ' .LoadConfig::getParam($rH_cfg, 'rnaQc','topTranscript');
    $command .= ' -s ' .$inputFile;
    $command .= ' -t ' .LoadConfig::getParam($rH_cfg, 'rnaQc','referenceGtf');
    $command .= ' -r ' .LoadConfig::getParam($rH_cfg, 'rnaQc','referenceFasta');
    $command .= ' -o ' .$outputFolder ;
    $command .= ' -BWArRNA ' .LoadConfig::getParam($rH_cfg, 'rnaQc','ribosomalFasta') .' &&';
    $command .= ' zip -r ' .$outputFolder .'.zip ' .$outputFolder;
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
  }
    
  return $ro_job;
}

sub saturation {
	my $rH_cfg      = shift;
	my $countFile   = shift;
	my $geneSizeFile     = shift;
	my $rpkmDir = shift;
	my $saturationDir = shift;


	my $up2date = PipelineUtils::testInputOutputs([$countFile], [$saturationDir .'.zip']);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
  	my $command;
		$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'saturation' , 'moduleVersion.cranR') .' ' . LoadConfig::getParam($rH_cfg, 'saturation' , 'moduleVersion.tools') . ' &&';
		$command .= ' Rscript \$R_TOOLS/rpkmSaturation.R ' .$countFile .' ' .$geneSizeFile .' ' .$rpkmDir .' ' .$saturationDir;
		$command .= ' ' .LoadConfig::getParam($rH_cfg, 'saturation' , 'threadNum');
		$command .= ' ' .LoadConfig::getParam($rH_cfg, 'saturation' , 'optionR').' &&';
		$command .= ' zip -r ' .$saturationDir .'.zip ' .$saturationDir;
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
	}
	return $ro_job;
}

sub fpkmCor {
	my $rH_cfg         = shift;
	my $paternFile     = shift;
	my $folderFile     = shift;
	my $outputBaseName = shift;
	
  my $up2date = PipelineUtils::testInputOutputs([$paternFile], [$outputBaseName]);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
  	my $command;
	  $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.cranR') .' ' . LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.tools') . ' ;';
  	$command .= ' Rscript \$R_TOOLS/fpkmStats.R ' .$paternFile .' ' .$folderFile .' ' .$outputBaseName;
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
	}
	return $ro_job;
}

sub readStats {
	my $rH_cfg = shift;
	my $inputFile = shift;
	my $outputFile = shift;
	my $sampleName = shift;
	my $fileType = shift;
	
  my $up2date = PipelineUtils::testInputOutputs([$inputFile], [$outputFile]);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
  	my $command;
		if ((lc $fileType) eq "trim") {
			$command .= 'grep \"Input Read\" ' .$inputFile;
			$command .= ' | awk -F\":\" -v na='.$sampleName .' \' BEGIN { OFS=\"\t\" } { split(\$2,a,\" \"); split(\$3,b,\" \"); print na,a[1],b[1]} \''; 
			$command .= ' > ' .$outputFile;
		}
		elsif  ((lc $fileType) eq "bam") {
			##decrepited
			my $unmapOption = '-F260';
			$command .= SAMtools::viewFilter($rH_cfg, $inputFile, $unmapOption, undef);
			$command .= ' | awk \' { print \$1} \'';
			#---- flefebvr Tue 16 Apr 13:47:53 2013 
			#$command .= ' | sort -u | wc -l  > ' .$outputFile;
			$command .= ' | sort -u -T '.LoadConfig::getParam($rH_cfg, 'metrics' , 'tmpDir');
			$command .= ' | wc -l  > ' .$outputFile;
			#----
		}
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
	}

	return $ro_job;
}


##Decrepated everthing should be in the report
sub mergeIndvidualReadStats{
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rawFile  = shift;
	my $filterFile = shift;
	my $alignFile = shift;
	my $outputFile =shift;
	
  my $up2date = PipelineUtils::testInputOutputs([$alignFile], [$outputFile]);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
	  my $command;
		$command .= 'echo \"' .$sampleName .'\"' ;
		$command .= ' | cat - ' .$rawFile .' ' .$filterFile .' ' .$alignFile;
		$command .= ' | tr \'\n\' \',\' ';
		$command .= ' | awk \' BEGIN {FS=\",\"} { print \$1 \",\" \$2 \",\" \$3 \",\" \$4}\''; 
		$command .= ' > ' .$outputFile .' &&'; 
		$command .= ' rm  ' .$rawFile .' ' .$filterFile .' ' .$alignFile;
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
	}

	return $ro_job;
}

sub mergeReadStats{
	my $rH_cfg         = shift;
	my $paternFile     = shift;
	my $folderFile     = shift;
	my $outputFile     = shift;

  my $up2date = PipelineUtils::testInputOutputs([$paternFile], [$outputFile]);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
	my $command;
  	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.cranR') .' ' . LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.tools') . ' ;';
	  $command .= ' Rscript \$R_TOOLS/mergeReadStat.R ' .$paternFile .' ' .$folderFile .' ' .$outputFile .' &&';
  	$command .= ' rm ' .$folderFile .'/*' .$paternFile;
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
	}
	return $ro_job;
}


sub mergeTrimmomaticStats{
	my $rH_cfg         = shift;
	my $libraryType    = shift;
	my $paternFile     = shift;
	my $folderFile     = shift;
	my $outputFile     = shift;

	if (!defined($libraryType) || $libraryType eq "") {
		$libraryType= 'unknown';
	}

  my $up2date = PipelineUtils::testInputOutputs([$paternFile], [$outputFile]);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
  	my $command;
	  $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.cranR') .' ' . LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.tools') . ' ;';
  	$command .= ' Rscript \$R_TOOLS/mergeTrimmomaticStat.R ' .$paternFile .' ' .$folderFile .' ' .$outputFile .' ' .$libraryType;
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
  }
	
	return $ro_job;
}

sub mergePrintReadStats{
	my $rH_cfg 							= shift;
	my $sampleName         = shift;
	my $trimOutputFile  		= shift;
	my $flagStatsFile 			= shift;
	my $outputFile  				= shift;
	my $command							= undef;

  my $up2date = PipelineUtils::testInputOutputs([$flagStatsFile], [$outputFile]);
  my $ro_job = new Job(!defined($up2date));

  if (!$ro_job->isUp2Date()) {
		$command .= 'module load ' . LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.tools') . ' ;';
		$command .= ' perl -MReadMetrics -e \' ReadMetrics::mergeStats(\"'.$sampleName.'\",';
		$command .= ' \"'. $outputFile .'\", ReadMetrics::parseTrimOutput(\"'.$sampleName.'\",';
		$command .= ' \"'. $trimOutputFile .'\"), ReadMetrics::parseFlagstats(\"'.$sampleName.'\",';
		$command .= ' \"'. $flagStatsFile .'\"));\'';
    $command .= ' ' . $up2date;

    $ro_job->addCommand($command);
	}
	return $ro_job;
}

1;
