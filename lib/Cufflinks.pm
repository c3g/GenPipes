#!/usr/env/perl

=head1 NAME

I<Cufflinks>

=head1 SYNOPSIS

Cufflinks-> fpkm()

=head1 DESCRIPTION

B<Cufflinks> is a library to do: Transcript assembly, differential expression, and differential regulation for RNA-Seq

Input = file_name

Output = array


=head1 AUTHOR
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package Cufflinks;

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
sub fpkm{
	my $rH_cfg        = shift;
	my $inputBAM      = shift;
	my $outputFolder     = shift;
	my $transcriptOption     = shift;

	my $latestFile = -M $inputBAM;
	my $outputIndexFile= $outputFolder. '/transcripts.gtf';

	my $command;
	# -M gives modified date relative to now. The bigger the older.
	if(!defined($latestFile) || !defined(-M $outputIndexFile) || $latestFile < -M $outputIndexFile) {
		$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'fpkm','cufflinksModule') .' ;';
		$command .= ' cufflinks -q';
		$command .= ' ' .$transcriptOption; 
		$command .= ' --max-bundle-frags ' .LoadConfig::getParam($rH_cfg, 'fpkm','cufflinksMaxFargs');
		$command .= ' --library-type ' .LoadConfig::getParam($rH_cfg, 'align', 'strandInfo');
		$command .= ' -p ' .LoadConfig::getParam($rH_cfg, 'fpkm','cufflinksThreads');
		$command .= ' -o ' .$outputFolder ;
		$command .= ' ' .$inputBAM;
	}
	
	return $command;
}

sub getDesign {
	my $rH_cfg        = shift;
	my $designFilePath = shitf;

	my %design;
	open(INFO, $designFilePath) or die "Can't find file $designFilePath\n";
	my @infos = <INFO>;
	close(INFO);

	my @splitA = split(/\t/, $infos[0]);
	my $numberDesigns = @splitA-1;
	for(my $i = 1; $i <= $numberDesigns; $i++) {

		my $designName = $splitA[$i];
		chomp($designName);
		my @group1;
		my @group2;

		for(my $j = 1; $j < @infos; $j++) {
		
			my @splitB = split(/\t/, $infos[$j]);
			my $sampleName = $splitB[0];
			chomp($sampleName);
			if($splitB[$i] == 1) {
				push(@group1,$sampleName);
			}
			elsif($splitB[$i] == 2) {
				 push(@group2,$sampleName); 
			}
			elsif($splitB[$i] == 0) {
			;	# do nothing
			}
			else {
				die "Wrong group assignment; check design file\n";
				
			}	

		}
		$design{$designName}=[@group1,@group2]
	}
	return \%design
}

1;
