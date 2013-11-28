#!/usr/env/perl

=head1 NAME

I<PacBioTools>

=head1 SYNOPSIS

PacBioTools->run()

=head1 DESCRIPTION

B<PacBioTools> This a library to analyze PacBio data using the SmrtAnalysis suite.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debug

=cut

package PacBioTools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Job;

# SUB
#-----------------------
sub getCutoff {
 	my $rH_cfg     = shift;
	my $infile     = shift;
	my $coverage   = shift;
	my $genomeSize = shift;
	my $xml        = shift;
	my $xmlOut     = shift;
	my $outfile    = shift;

  	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], [$outfile.".fasta"]);

	if (!$ro_job->isUp2Date()) {
		my $cmd;
	
		# Choose a subread length threshold such that subreads above the threshold provide about 20x coverage of the genome.
		$cmd = 'module load '.LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.mugqictools').' ;';
		$cmd .= ' pacBioGetCutoff.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --coverage ' . $coverage;
		$cmd .= ' --genomeSize ' . $genomeSize;
		$cmd .= ' --coverageCutoff';
		$cmd .= ' --coverageFraction ' . LoadConfig::getParam($rH_cfg, 'preassembly', 'coverageFraction');
		$cmd .= ' --xml ' .$xml;
		$cmd .= ' --xmlOut ' . $xmlOut;
		$cmd .= ' > '.$outfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub celeraConfig {
 	my $rH_cfg         = shift;
	my $merSize        = shift;
	my $infile         = shift;
	my $outfile        = shift;

  	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], [$outfile.".fasta"]);

	if (!$ro_job->isUp2Date()) {
		my $cmd;

		$cmd = 'module load '.LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.mugqictools').' ;';
		$cmd .= ' pacBioAssemblyCeleraConfig.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'num_threads');
		$cmd .= ' --minReadSize ' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'minReadSize');
		$cmd .= ' --overlapper ' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'overlapper');
		$cmd .= ' --merCompression ' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'merCompression');
		#$cmd .= ' --merSize ' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'merSize');
		$cmd .= ' --merSize ' . $merSize;
		$cmd .= ' --merylMemory ' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'merylMemory');
		$cmd .= ' > ' . $outfile;
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub assemblyStats{
 	my $rH_cfg                = shift;
	#my $indir                 = shift;
	my $filteredSubreadsTable = shift;
	my $filteredSummary       = shift;
	my $assemblyQc            = shift;
	my $sampleName            = shift;
	my $suffix                = shift;
	my $estimatedGenomeSize   = shift;
	my $smrtCells             = shift;
	my $outdir                = shift;

  	my $ro_job = new Job();
	$ro_job->testInputOutputs(
		[$filteredSubreadsTable, $filteredSubreadsTable, $filteredSummary],
		#[$indir."/filtering/data/filtered_summary.csv", $indir."/filtering/results/filterReports_filterStats.xml"], 
		[$outdir."/summaryTableAssembly.tsv", $outdir."/summaryTableReads.tsv"]
	);

	if (!$ro_job->isUp2Date()) {
		my $cmd;

		$cmd .= 'module load '.LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.R').' ;';
		$cmd .= 'module load '.LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.mugqictools').' ;';
		$cmd .= ' pacBioAssemblyStats.pl';
		#$cmd .= ' --indir ' . $indir;
		$cmd .= ' --filteredSubreadsTable ' . $filteredSubreadsTable;
		$cmd .= ' --filteredSummary ' . $filteredSummary;
		$cmd .= ' --assemblyQc ' . $assemblyQc;
		$cmd .= ' --sampleName ' . $sampleName;
		$cmd .= ' --suffix ' . $suffix;
		$cmd .= ' --estimatedGenomeSize ' . $estimatedGenomeSize;
		$cmd .= ' --smrtCells ' . $smrtCells;
		$cmd .= ' --outdir ' . $outdir;

		$ro_job->addCommand($cmd);

	}
	return $ro_job;
}

sub splitReads{
 	my $rH_cfg              = shift;
	my $subreads            = shift;
	my $cutoff              = shift; # a file containing the cutoff number is actually passed in arg here.
	my $shortReads          = shift;
	my $longReads           = shift;

  	my $ro_job = new Job();
	$ro_job->testInputOutputs([$subreads], [$shortReads, $longReads]);

	if (!$ro_job->isUp2Date()) {
		my $cmd;

		$cmd = 'module load '.LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.mugqictools').' ;';
		$cmd .= ' pacBioSplitReads.pl';
		$cmd .= ' --infile ' . $subreads;
		$cmd .= ' --cutoff `cat ' . $cutoff . '` ';
		$cmd .= ' --outfileShort ' . $shortReads;
		$cmd .= ' --outfileLong ' . $longReads;

		$ro_job->addCommand($cmd);

	}
	return $ro_job;
}

1;
