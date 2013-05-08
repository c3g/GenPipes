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
	
	if(!defined($transcriptOption)) {
		$transcriptOption = '';
	}
	my $command;
	# -M gives modified date relative to now. The bigger the older.
	if(!defined($latestFile) || !defined(-M $outputIndexFile) || $latestFile < -M $outputIndexFile) {
		$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'fpkm','moduleVersion.cufflinks') .' ;';
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
	my $rH_cfg  = shift;
	my $designFilePath = shift;

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
		$design{$designName}=[\@group1,\@group2];
	}
	return \%design;
}

sub cuffdiff {
	my $rH_cfg        = shift;
	my $rA_groupInputFiles        = shift;
	my $outputDir     = shift;
	my $referenceGtf     = shift;
	
	my $groupCmd;
	my $numberGroup = @{$rA_groupInputFiles};
	for (my $i=0 ; $i <  $numberGroup; $i++) {
		$groupCmd .= ' ' .$rA_groupInputFiles->[$i];
	}
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'cuffdiff','moduleVersion.cufflinks') .' ;';
	$command .= ' cuffdiff -p' .LoadConfig::getParam($rH_cfg, 'cuffdiff','numThreads');
	$command .= ' -o ' .$outputDir;
	$command .= ' ' .$referenceGtf;
	$command .= ' ' .$groupCmd;
	$command .= ' ' .LoadConfig::getParam($rH_cfg, 'cuffdiff','options');
	$command .= ' -b ' .LoadConfig::getParam($rH_cfg, 'cuffdiff','referenceFasta');

	return $command;
}

sub cuffmerge {
	my $rH_cfg        = shift;
	my $mergeListFile = shift;
	my $outputDir     = shift;
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'cuffmerge','moduleVersion.cufflinks') .' ;';
	$command .= ' cuffmerge -p' .LoadConfig::getParam($rH_cfg, 'cuffmerge','numThreads');
	$command .= ' -o ' .$outputDir;
	$command .= ' -g ' .LoadConfig::getParam($rH_cfg, 'cuffmerge','referenceGtf');
	$command .= ' -s ' .LoadConfig::getParam($rH_cfg, 'fpkm','referenceFasta');
	$command .= ' ' .$mergeListFile;

	return $command;
}


sub mergeGtfFormat {
	my $rH_cfg        = shift;
	my $inputFile     = shift;
	my $outputFile    = shift;
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'cuffmerge','moduleVersion.tools') .' ;';
	$command .= ' perl $PERL_TOOLS/formatGtfCufflinks.pl' .' '.$inputFile .' ' .$outputFile ;
	
	return $command;
}

sub mergeCuffdiffRes {
	my $rH_cfg        = shift;
	my $designFile    = shift;
	my $outputDir     = shift;

	### TO DO : re-write mergecuffdiff_known.R and mergecuffdiff_denovo.R to be more portable
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'cuffmerge','moduleVersion.tools') .' ' .LoadConfig::getParam($rH_cfg, 'cuffdiff','moduleVersion.cranR') .' ;';
	$command .= ' Rscript $R_TOOLS/mergecuffdiff_known.R ' .$outputDir .$designFile .' &&';
	$command .= ' Rscript $R_TOOLS/mergecuffdiff_denovo.R ' .$outputDir .$designFile ;

	return $command;
}

sub filterResults {
	my $rH_cfg        = shift;
	my $outputDir     = shift;
	
	### TO DO : make it more portable when the mergeCuffdiffRes R script will be re-write
	my $command;
	$command .= 'for i in \`ls ' .$outputDir .'/*/isoform_exp.diff.with.fpkm.csv\` ;' ;
	$command .= ' do head -1 \$i > ' .$outputDir .'/tmp ;' ;
	$command .= ' sed 1d \$i | grep -v \"NOTEST\" | grep -v \"FAIL\" | sort -k 12 -g >> ' .$outputDir .'/tmp ;' ;
	$command .= ' mv ' .$outputDir .'/tmp \$i ;' ;
	$command .= ' done' ;

	return $command;
}

1;
