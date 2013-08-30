#!/usr/env/perl

=head1 NAME

I<Breakdancer>

=head1 SYNOPSIS

Breakdancer

=head1 DESCRIPTION

B<Breakdancer> is a library to analyse SV events in a genome

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package Breakdancer;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use PipelineUtils;

# SUB
#-----------------------
sub bam2cfg {
    my $rH_cfg        = shift;
    my $sampleName    = shift;
    my $sampleBAM     = shift;
    my $output        = shift;
    my $stdDevCutoff  = shift;

    if(!defined($stdDevCutoff)) {
      $stdDevCutoff = 3;
    }
  
    my $up2date = PipelineUtils::testInputOutputs([$sampleBAM], [$output]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'bam2cfg', 'moduleVersion.breakdancer').' ;';
      $command .= ' module load '.LoadConfig::getParam($rH_cfg, 'bam2cfg', 'moduleVersion.samtools').' ;';
      $command .= ' bam2cfg.pl -g -h ';
      $command .= ' -c '.$stdDevCutoff;
      $command .= ' '.$sampleBAM;
      $command .= ' > '.$output;
      $command .= ' ' . $up2date;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub pairedBRDITX {
    my $rH_cfg          = shift;
    my $sampleName      = shift;
    my $inputCFG        = shift;
    my $outputPrefix    = shift;

    my $up2date = PipelineUtils::testInputOutputs([$inputCFG], [$outputPrefix.'.bed', $outputPrefix.'.ctx']);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pairedBRDITX', 'moduleVersion.breakdancer').' ;';
      $command .= ' breakdancer_max '.LoadConfig::getParam($rH_cfg, 'pairedBRDITX', 'brdParameters');
      $command .= ' -g '.$outputPrefix.'.bed';
      $command .= ' -d '.$outputPrefix.'.ctx';
      $command .= ' '.$inputCFG;
      $command .= ' > '.$outputPrefix.'.ctx';
      $command .= ' ' . $up2date;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub pairedBRD {
    my $rH_cfg          = shift;
    my $sampleName      = shift;
    my $chr             = shift;
    my $inputCFG        = shift;
    my $outputPrefix    = shift;

    my $up2date = PipelineUtils::testInputOutputs([$inputCFG], [$outputPrefix.'.bed', $outputPrefix.'.ctx']);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pairedBRD', 'moduleVersion.breakdancer').' ;';
      $command .= ' breakdancer_max '.LoadConfig::getParam($rH_cfg, 'pairedBRD', 'brdParameters');
      $command .= ' -o '.$chr;
      $command .= ' -g '.$outputPrefix.'.bed';
      $command .= ' -d '.$outputPrefix.'.ctx';
      $command .= ' '.$inputCFG;
      $command .= ' > '.$outputPrefix.'.ctx';
      $command .= ' ' . $up2date;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub mergeCTX {
    my $rH_cfg          = shift;
    my $outputPrefix    = shift;

    my $up2date = PipelineUtils::testInputOutputs(undef, [$outputPrefix .'.ctx']);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'rm ' .$outputPrefix .'.ctx && ' ;
      $command .= 'touch ' .$outputPrefix .'.ctx && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*.ctx ; do cat \$i >> '  .$outputPrefix .'.ctx' ;
      $command .= ' ' . $up2date;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}
 
1;
