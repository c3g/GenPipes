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
    my $tumorPileup  = shift;
    my $normalPileup = shift;
    my $sampleConfig  = shift;
    my $output        = shift;

    my $outDate = -M $sampleConfig;
    my $inDate = -M $normalPileup;
    my $inDate2 = -M $tumorPileup;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    if(!defined($outDate) || !defined($inDate) || $inDate < $outDate || !defined($inDate2) || $inDate2 < $outDate ) {
        $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'moduleVersion.controlFreec') .' && ';
        $command .= 'sed \"s|TUMOR_PILEUP|' .$tumorPileup .'|g\" ' .'\${FREEC_HOME}/' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'referenceConfigFile') .' > ' .$sampleConfig .' && ';
        $command .= 'sed \"s|NORMAL_PILEUP|' .$normalPileup .'|g\" -i ' .$sampleConfig .' && ';
        $command .= 'sed \"s|OUTPUT_DIR|' .$output .'|g\"  -i ' .$sampleConfig .' && ';
        $command .= 'sed \"s|FORMAT_TYPE|' .LoadConfig::getParam($rH_cfg, 'ControlFreec', 'inputType') .'|g\"  -i ' .$sampleConfig .' && ';
        $command .= 'freec -conf ' .$sampleConfig ;
    }
    
    return $command;
}

    
1;
