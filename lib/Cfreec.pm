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
use File::Basename;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub pairedFreec {
    my $rH_cfg        = shift;
    my $tumorBam  = shift;
    my $normalBam = shift;
    my $sampleConfig  = shift;
    my $output        = shift;

    my $outDate = -M $sampleConfig;
    my $inDate = -M $normalBam;
    my $inDate2 = -M $tumorBam;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    if(!defined($outDate) || !defined($inDate) || $inDate < $outDate || !defined($inDate2) || $inDate2 < $outDate ) {
        $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'moduleVersion.controlFreec') .' ' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'moduleVersion.samtools') .' && ';
        $command .= 'sed \"s|TUMOR_PILEUP|' .$tumorBam .'|g\" ' .'\${FREEC_HOME}/' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'referenceConfigFile') .' > ' .$sampleConfig .' && ';
        $command .= 'sed \"s|NORMAL_PILEUP|' .$normalBam .'|g\" -i ' .$sampleConfig .' && ';
        $command .= 'sed \"s|OUTPUT_DIR|' .$output .'|g\"  -i ' .$sampleConfig .' && ';
        $command .= 'freec -conf ' .$sampleConfig ;
    }
    
    return $command;
}

    
1;
