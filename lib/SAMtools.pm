#!/usr/env/perl

=head1 NAME

I<SAMtools>

=head1 SYNOPSIS

SAMtools->mpileup()

=head1 DESCRIPTION

B<SAMtools> is a library to manipulate and compute stats off of BAMs and to call variants

Input = file_name

Output = array

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package SAMtools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub mpileupPaired {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $normalBam   = shift;
  my $tumorBam    = shift;
  my $seqName     = shift;
  my $outputDir   = shift;

  my @bams = ($normalBam, $tumorBam);

  return mpileupBuilder($rH_cfg, $seqName, $outputDir, $sampleName, \@bams, 1);
}

sub mpileup {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rA_bams     = shift;
  my $seqName     = shift;
  my $outputDir   = shift;

  return mpileupBuilder($rH_cfg, $seqName, $outputDir, $sampleName, $rA_bams);
}

sub mpileupBuilder {
  my $rH_cfg      = shift;
  my $seqName     = shift;
  my $outputDir   = shift;
  my $sampleName  = shift;
  my $rA_bams     = shift;
  my $isPaired    = shift;

  if(defined($isPaired) && $isPaired == 1 && @{$rA_bams} != 2) {
    die("Paired was asked but there aren't 2 bams given\n");
  }

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $outputBCF = $outputDir.$sampleName.'.'.$seqName.'.bcf'; 

  my $command;
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'mpileup', 'moduleVersion.samtools').' ;';
  $command .= ' samtools mpileup';
  $command .= ' '.LoadConfig::getParam($rH_cfg, 'mpileup', 'mpileupExtraFlags');
  $command .= ' -f '.$refGenome;
  $command .= ' -r '.$seqName;
  for my $bamFiles (@{$rA_bams}) {
    $command .= ' '.$bamFiles;
  }
  if(defined($isPaired) && $isPaired == 1) {
    $command .= ' | bcftools view -T pair -bvcg - > '.$outputBCF;
  }
  else {
    $command .= ' | bcftools view -bvcg - > '.$outputBCF;
  }

  return $command;
}

sub mergeFilterBCF {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $bcfDir      = shift;
  my $outputDir   = shift; 
  my $rA_seqNames = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $outputBCF = $outputDir.$sampleName.'.merged.bcf'; 
  my $outputVCF = $outputDir.$sampleName.'.merged.flt.vcf'; 

  my $bcfInputs = "";
  for my $seqName (@$rA_seqNames) {
    my $bcfFile = $bcfDir.$sampleName.'.'.$seqName.'.bcf'; 

    $bcfInputs .= $bcfFile.' ';
  }
  
  my $command;
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'mergeFilterBCF', 'moduleVersion.samtools').' ;';
  $command .= ' bcftools cat';
  $command .= ' '.$bcfInputs;
  $command .= ' > '.$outputBCF;
  $command .= ' ; bcftools view '.$outputBCF;
  $command .= ' | vcfutils.pl varFilter';
  $command .= ' '.LoadConfig::getParam($rH_cfg, 'mergeFilterBCF', 'varfilterExtraFlags');
  $command .= ' > '.$outputVCF;

  return $command;
}

sub flagstat {
  my $rH_cfg     = shift;
  my $bamFile    = shift;
  my $output     = shift;

  my $command;
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'flagstat', 'moduleVersion.samtools').' ;';
  $command .= ' samtools flagstat';
  $command .= ' '.$bamFile;
  $command .= ' > '.$output;
}

sub idxstats {
  my $rH_cfg     = shift;
  my $bamFile    = shift;
  my $output     = shift;

  my $command;
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'idxstats', 'moduleVersion.samtools').' ;';
  $command .= ' samtools idxstats';
  $command .= ' '.$bamFile;
  $command .= ' > '.$output;
}

sub rawmpileup {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $bamFile     = shift;
  my $seqName     = shift;
  my $output      = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');

  my $command;
  $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'rawmpileup', 'moduleVersion.samtools').' ;';
  $command .= ' samtools mpileup';
  $command .= ' '.LoadConfig::getParam($rH_cfg, 'rawmpileup', 'mpileupExtraFlags');
  $command .= ' -f '.$refGenome;
  $command .= ' -r '.$seqName;
  $command .= ' '.$bamFile;
  $command .= ' | gzip -1 -c > '.$output;
}

sub sort {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $option      = shift;
  my $bamFile     = shift;
  my $output      = shift;

 
  my $command;
  $command .= 'module load ' .LoadConfig::getParam($rH_cfg,'default','moduleVersion.samtools') .' ;';
  $command .= ' samtools sort';
  $command .= ' '.$option;
  $command .= ' '.$bamFile;
  $command .= ' '.$output;

  return $command;
}

sub viewFilter {
	my $rH_cfg      = shift;
	my $bamFile     = shift;
	my $option      = shift;
	my $output      = shift;

	my $returnOutput = '';
	if (defined($output)) {
		$returnOutput = ' > ' .$output;
	}
	if (!(defined($option))) {
		$option = '';
	}
	
	
	my $command;
	$command .= 'module load ' .LoadConfig::getParam($rH_cfg,'default','moduleVersion.samtools') .' ;';
	$command .= ' samtools view';
	$command .= ' ' .$option;
	$command .= ' ' .$bamFile;
	$command .= $returnOutput;

	return $command;
}

1;
