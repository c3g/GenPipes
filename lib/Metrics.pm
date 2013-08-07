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
use SAMtools;

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
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'rnaQc','moduleVersion.bwa') ;
    $command .= ' ' .LoadConfig::getParam($rH_cfg, 'rnaQc','moduleVersion.rnaseqc') .' &&';
    $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'rnaQc', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'rnaQc', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'rnaQc', 'metricsRam').' -jar \${RNASEQC_JAR}';
    $command .= ' -n ' .LoadConfig::getParam($rH_cfg, 'rnaQc','topTranscript');
    $command .= ' -s ' .$inputFile;
    $command .= ' -t ' .LoadConfig::getParam($rH_cfg, 'rnaQc','referenceGtf');
    $command .= ' -r ' .LoadConfig::getParam($rH_cfg, 'rnaQc','referenceFasta');
    $command .= ' -o ' .$outputFolder ;
    $command .= ' -BWArRNA ' .LoadConfig::getParam($rH_cfg, 'rnaQc','ribosomalFasta');
  }
    
  return $command;
}

sub saturation {
	my $rH_cfg      = shift;
	my $countFile   = shift;
	my $geneSizeFile     = shift;
	my $rpkmDir = shift;
	my $saturationDir = shift;


	

	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'saturation' , 'moduleVersion.cranR') .' ' . LoadConfig::getParam($rH_cfg, 'saturation' , 'moduleVersion.tools') . ' &&';
	$command .= ' Rscript \$R_TOOLS/rpkmSaturation.R ' .$countFile .' ' .$geneSizeFile .' ' .$rpkmDir .' ' .$saturationDir;
	$command .= ' ' .LoadConfig::getParam($rH_cfg, 'saturation' , 'threadNum');
	$command .= ' ' .LoadConfig::getParam($rH_cfg, 'saturation' , 'optionR');

	return $command;
}

sub fpkmCor {
	my $rH_cfg         = shift;
	my $paternFile     = shift;
	my $folderFile     = shift;
	my $outputBaseName = shift;
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.cranR') .' ' . LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.tools') . ' ;';
	$command .= ' Rscript \$R_TOOLS/fpkmStats.R ' .$paternFile .' ' .$folderFile .' ' .$outputBaseName;
	
	return $command;
}

sub readStats {
	my $rH_cfg = shift;
	my $inputFile = shift;
	my $outputFile = shift;
	my $sampleName = shift;
	my $fileType = shift;
	
	my $latestInputFile = -M $inputFile;
	my $latestOutputFile = -M $outputFile;
	
	my $command;
	if(!defined($latestInputFile) || !defined($latestOutputFile) || $latestInputFile <  $latestOutputFile) {
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
	}

	return $command;
}


##Decrepated everthing should be in the report
sub mergeIndvidualReadStats{
	my $rH_cfg = shift;
	my $sampleName = shift;
	my $rawFile  = shift;
	my $filterFile = shift;
	my $alignFile = shift;
	my $outputFile =shift;
	
	my $latestInputFile = -M $alignFile;
	my $latestOutputFile = -M $outputFile;
	
	my $command;
	if(!defined($latestInputFile) || !defined($latestOutputFile) || $latestInputFile <  $latestOutputFile) {
		$command .= 'echo \"' .$sampleName .'\"' ;
		$command .= ' | cat - ' .$rawFile .' ' .$filterFile .' ' .$alignFile;
		$command .= ' | tr \'\n\' \',\' ';
		$command .= ' | awk \' BEGIN {FS=\",\"} { print \$1 \",\" \$2 \",\" \$3 \",\" \$4}\''; 
		$command .= ' > ' .$outputFile .' &&'; 
		$command .= ' rm  ' .$rawFile .' ' .$filterFile .' ' .$alignFile;
	}

	return $command;
}

sub mergeReadStats{
	my $rH_cfg         = shift;
	my $paternFile     = shift;
	my $folderFile     = shift;
	my $outputFile     = shift;
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.cranR') .' ' . LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.tools') . ' ;';
	$command .= ' Rscript \$R_TOOLS/mergeReadStat.R ' .$paternFile .' ' .$folderFile .' ' .$outputFile .' &&';
	$command .= ' rm ' .$folderFile .'/*' .$paternFile;
	
	return $command;
}


sub mergePrintReadStats{
	my $rH_cfg 							= shift;
	my $sampleName         = shift;
	my $trimOutputFile  		= shift;
	my $flagStatsFile 			= shift;
	my $outputFile  				= shift;
	my $command							= undef;

	my $latestInputFile = -M $flagStatsFile;
	my $latestOutputFile = -M $outputFile;
	
	if(!defined($latestInputFile) || !defined($latestOutputFile) || $latestInputFile <  $latestOutputFile) {
		$command .= 'module load ' . LoadConfig::getParam($rH_cfg, 'metrics' , 'moduleVersion.tools') . ' ;';
		$command .= ' perl -MReadMetrics -e \' ReadMetrics::mergeStats(\"'.$sampleName.'\",';
		$command .= ' \"'. $outputFile .'\", ReadMetrics::parseTrimOutput(\"'.$sampleName.'\",';
		$command .= ' \"'. $trimOutputFile .'\"), ReadMetrics::parseFlagstats(\"'.$sampleName.'\",';
		$command .= ' \"'. $flagStatsFile .'\"));\'';
	}
	return $command;
}

1;
