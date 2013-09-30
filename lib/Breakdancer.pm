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
  
    my $ro_job = new Job();
    $ro_job->testInputOutputs([$sampleBAM], [$output]);

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'bam2cfg', 'moduleVersion.breakdancer').' ;';
      $command .= ' module load '.LoadConfig::getParam($rH_cfg, 'bam2cfg', 'moduleVersion.samtools').' ;';
      $command .= ' bam2cfg.pl -g -h ';
      $command .= ' -c '.$stdDevCutoff;
      $command .= ' '.$sampleBAM;
      $command .= ' > '.$output;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub pairedBRDITX {
    my $rH_cfg          = shift;
    my $sampleName      = shift;
    my $inputCFG        = shift;
    my $outputPrefix    = shift;

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$inputCFG], [$outputPrefix.'.bed', $outputPrefix.'.ctx']);

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pairedBRDITX', 'moduleVersion.breakdancer').' ;';
      $command .= ' breakdancer_max '.LoadConfig::getParam($rH_cfg, 'pairedBRDITX', 'brdParameters');
      $command .= ' -g '.$outputPrefix.'.bed';
      $command .= ' -d '.$outputPrefix.'.ctx';
      $command .= ' '.$inputCFG;
      $command .= ' > '.$outputPrefix.'.ctx';

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

    my $ro_job = new Job();
    $ro_job->testInputOutputs([$inputCFG], [$outputPrefix.'.bed', $outputPrefix.'.ctx']);

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pairedBRD', 'moduleVersion.breakdancer').' ;';
      $command .= ' breakdancer_max '.LoadConfig::getParam($rH_cfg, 'pairedBRD', 'brdParameters');
      $command .= ' -o '.$chr;
      $command .= ' -g '.$outputPrefix.'.bed';
      $command .= ' -d '.$outputPrefix.'.ctx';
      $command .= ' '.$inputCFG;
      $command .= ' > '.$outputPrefix.'.ctx';

      $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub mergeCTX {
    my $rH_cfg          = shift;
    my $outputPrefix    = shift;

    my $ro_job = new Job();
    $ro_job->testInputOutputs(undef, [$outputPrefix .'.ctx']);

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'rm ' .$outputPrefix .'.ctx && ' ;
      $command .= 'touch ' .$outputPrefix .'.ctx && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*.ctx ; do cat \$i >> '  .$outputPrefix .'.ctx' ;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}
 
1;
