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

package Picard;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub runPairedDNAC {
    my $rH_cfg        = shift;
    my $sampleName    = shift;
    my $inputBins     = shift;
    my $outputPrefix  = shift;
    my $window        = shift;

    my $outDate = -M $outputPrefix.'.txt';
    my $inDate = -M $inputBins;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'runPairedDNAC', 'moduleVersion.cranR').' ;';
        $command .= ' module load '.LoadConfig::getParam($rH_cfg, 'runPairedDNAC', 'moduleVersion.svtools').' ;';
        $command .= ' Rscript \${SVTOOLS_HOME}/Cancer/RunDNAC.6.0.R';
        $command .= ' -f '.$inputBins;
        $command .= ' -b '.$window;
        $command .= ' -o '.$outputPrefix;
    }
    return $command;
}

sub filterDNAC {
    my $rH_cfg          = shift;
    my $sampleName      = shift;
    my $inputDNACCalls  = shift;
    my $outputPrefix    = shift;
    my $cnvProx           = shift;

    my $outDate = -M $outputPrefix.'.filteredSV.txt';
    my $inDate = -M $inputDNACCalls;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'runPairedDNAC', 'moduleVersion.svtools').' ;';
        $command .= ' \${SVTOOLS_HOME}/Cancer/filterOutDNAC.sh';
        $command .= ' '.$inputDNACCalls;
        $command .= ' '.$outputPrefix.'.filteredSV.txt';
        $command .= ' '.$sampleName;
        $command .= ' 5';
        $command .= ' && ';
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'runPairedDNAC', 'moduleVersion.svtools').' ;';
        $command .= ' \${SVTOOLS_HOME}/Cancer/filterBedResults.sh';
        $command .= ' '.$outputPrefix.'.filteredSV.txt';
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'filterDNAC', 'referenceMappabilityBed');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'filterDNAC', 'referenceGeneCoordinates');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'filterDNAC', 'referenceDGVCoordinates');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'filterDNAC', 'referenceMicrosatellitesCoordinates');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'filterDNAC', 'referenceRepeatMaskerCoordinates');
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'filterDNAC', 'referenceGenomeLengths');
        $command .= ' '.$outputPrefix;
        $command .= ' '.$outputPrefix.'.tmp';
        $command .= ' '.$cnvProx;
    }
    return $command;
}
1;
