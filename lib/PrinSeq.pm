#!/usr/env/perl

=head1 NAME

I<Picard>

=head1 SYNOPSIS

SmrtAnalysis->run()

=head1 DESCRIPTION

B<Picard> This a library to analyze PAcBio data using the SmrtAnalysis suite.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package PrinSeq;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub fastqToFasta {
 	my $rH_cfg     = shift;
	my $infile     = shift;
	my $outfile    = shift;

  	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], [$outfile.".fasta"]);

	if (!$ro_job->isUp2Date()) {
		my $cmd;
		# Fastq to Fasta 
		$cmd = 'module load '.LoadConfig::getParam($rH_cfg, 'prinseq', 'moduleVersion.prinseq').' ;';
		$cmd .= ' prinseq-lite.pl';
		$cmd .= ' -verbose';
		$cmd .= ' -fastq ' . $infile;
		$cmd .= ' -out_format 1';
		$cmd .= ' -out_good ' . $outfile;
		$cmd .= ' > /dev/null 2>&1';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}


1;
