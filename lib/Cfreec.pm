#!/usr/env/perl

=head1 NAME

I<Pindel>

=head1 SYNOPSIS

Pindel

=head1 DESCRIPTION

B<Cfreec> is a library to analyse CNV events in a genome

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package Cfreec;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use File::Basename;
use PipelineUtils;

# SUB
#-----------------------
sub pairedFreec {
    my $rH_cfg        = shift;
    my $tumorPileup  = shift;
    my $normalPileup = shift;
    my $sampleConfig  = shift;
    my $output        = shift;

    my $up2date = PipelineUtils::testInputOutputs([$normalPileup, $tumorPileup], [$sampleConfig]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'moduleVersion.controlFreec') .' && ';
        $command .= 'sed \"s|TUMOR_PILEUP|' .$tumorPileup .'|g\" ' .'\${FREEC_HOME}/' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'referenceConfigFile') .' > ' .$sampleConfig .' && ';
        $command .= 'sed \"s|NORMAL_PILEUP|' .$normalPileup .'|g\" -i ' .$sampleConfig .' && ';
        $command .= 'sed \"s|OUTPUT_DIR|' .$output .'|g\"  -i ' .$sampleConfig .' && ';
        $command .= 'sed \"s|FORMAT_TYPE|' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'inputType') .'|g\"  -i ' .$sampleConfig .' && ';
        $command .= 'freec -conf ' .$sampleConfig ;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }
    
    return $ro_job;
}

    
1;
