#!/usr/bin/env perl

=head1 NAME

I<SVTools>

=head1 SYNOPSIS

Picard->merge()

=head1 DESCRIPTION

B<SVTools> is a library to analyse BAMs for Structural Variants

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package VCFtools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin";

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------
sub annotateMappability {
  my $rH_cfg    = shift;
  my $inputVCF  = shift;
  my $outputVCF = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputVCF], [$outputVCF]);

  if (!$ro_job->isUp2Date()) {
    my $command;

    $command .= LoadConfig::moduleLoad($rH_cfg, [
      ['annotateMappability', 'moduleVersion.vcftools'],
      ['annotateMappability', 'moduleVersion.tabix']
    ]) . ' &&';
    $command .= ' vcf-annotate -d \"key=INFO,ID=MIL,Number=1,Type=String,Description=' . "'" . 'Mappability annotation. 300IS 40SD 1SHI. HC = to high coverage (>400), LC = to high coverage (<50), MQ = to low mean mapQ (<20), ND = no data at the position' . "'" . '\"';
    $command .= ' -c CHROM,FROM,TO,INFO/MIL ';
    $command .= ' -a ' . LoadConfig::getParam($rH_cfg, 'annotateMappability', 'referenceMappabilityBedIndexed', 1, 'filepath');
    $command .= ' ' . $inputVCF;
    $command .= ' > ' . $outputVCF;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub indexVCF {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $inputVCF    = shift;
  my $outputVCF   = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputVCF], [$outputVCF]);

  if (!$ro_job->isUp2Date()) {
    my $command;

    $command .= LoadConfig::moduleLoad($rH_cfg, [['indexVCF', 'moduleVersion.tabix']]) . ' &&';
    $command .= ' bgzip -c';
    $command .= ' ' . $inputVCF;
    $command .= ' > ' . $outputVCF;
    $command .= ' && ';
    $command .= ' tabix -p vcf -f';
    $command .= ' ' . $outputVCF;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

sub mergeVCF {
  my $rH_cfg    = shift;
  my $rA_vcfs   = shift;
  my $outputVCF = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs($rA_vcfs, [$outputVCF]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['mergeVCF', 'moduleVersion.vcftools']]) . ' &&';
    $command .= ' vcf-concat';
    $command .= ' ' . join(' ', @{$rA_vcfs});
    $command .= ' > ' . $outputVCF;

    $ro_job->addCommand($command);
  }
  return $ro_job;
}

1;
