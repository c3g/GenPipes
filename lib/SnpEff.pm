#!/usr/env/perl

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

package SnpEff;

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
sub annotateDbSnp {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $inputVCF    = shift;
    my $outputVCF   = shift;

    my $up2date = PipelineUtils::testInputOutputs([$inputVCF],[$outputVCF]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {  
        my $command;
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'moduleVersion.snpeff').' ;';
        $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'siftRam');
        $command .= ' -jar \${SNPEFF_HOME}/SnpSift.jar annotate';
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'annotateDbSnp', 'dbSnp');
        $command .= ' '.$inputVCF;
        $command .= ' > '.$outputVCF;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub computeEffects {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $inputVCF    = shift;
    my $outputVCF   = shift;
  
    my $up2date = PipelineUtils::testInputOutputs([$inputVCF],[$outputVCF]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'computeEffects', 'moduleVersion.snpeff').' ;';
        $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'computeEffects', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'computeEffects', 'extraJavaFlags'). ' -Xmx'.LoadConfig::getParam($rH_cfg, 'computeEffects', 'snpeffRam');
        $command .= ' -jar \${SNPEFF_HOME}/snpEff.jar eff';
        $command .= ' -c \${SNPEFF_HOME}/snpEff.config';
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'computeEffects', 'snpeffParams');
        $command .= ' -o vcf';
        $command .= ' -i vcf';
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'computeEffects', 'referenceSnpEffGenome');
        $command .= ' '.$inputVCF;
        $command .= ' > '.$outputVCF;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }
    return $ro_job;
}

sub annotateDbNSFP {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $inputVCF    = shift;
    my $outputVCF   = shift;
  
    my $up2date = PipelineUtils::testInputOutputs([$inputVCF],[$outputVCF]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'moduleVersion.snpeff').' ;';
        $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'siftRam');
        $command .= ' -jar \${SNPEFF_HOME}/SnpSift.jar dbnsfp';
        $command .= ' -v '.LoadConfig::getParam($rH_cfg, 'annotateDbNSFP', 'dbNSFP');
        $command .= ' '.$inputVCF;
        $command .= ' > '.$outputVCF;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }
    return $ro_job;
}

1;
