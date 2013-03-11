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

package SVtools;

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
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'runPairedDNAC', 'moduleVersion.cranR').' ;';
        $command .= ' module load '.LoadConfig::getParam($rH_cfg, 'runPairedDNAC', 'moduleVersion.svtools').' ;';
        $command .= ' Rscript \${SVTOOLS_HOME}/Cancer/RunDNAC.6.0.R';
        $command .= ' -f '.$inputBins;
        $command .= ' -b '.$window;
        $command .= ' -o '.$outputPrefix;
    #}
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
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'filterSV', 'moduleVersion.svtools').' ;';
        $command .= ' \${SVTOOLS_HOME}/Cancer/filterOutDNAC.sh';
        $command .= ' '.$inputDNACCalls;
        $command .= ' '.$outputPrefix.'.txt';
        $command .= ' '.$sampleName;
        $command .= ' '.LoadConfig::getParam($rH_cfg, 'filterSV', 'minBinCNV');
        $command .= ' && ';
        $command .= filterResult($rH_cfg,  'filterSV', $outputPrefix, $cnvProx) ;
        $command .= ' && ';
        $command .= generateBedResult($rH_cfg, $outputPrefix) ;
    #}
    return $command;
}


sub filterBrD {
    my $rH_cfg          = shift;
    my $sampleName      = shift;
    my $inputBrDCalls  = shift;
    my $outputPrefix    = shift;
    my $normalFile    = shift;
    my $tumorFile    = shift;

    my $outDate = -M $outputPrefix.'.filteredSV.txt';
    my $inDate = -M $inputDNACCalls;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'filterSVC', 'moduleVersion.svtools').' ;';
        $command .= ' \${SVTOOLS_HOME}/Cancer/filterOutBrD.py';
        $command .= ' -f ' .$inputBrDCalls;
        $command .= ' -o ' .$outputPrefix.'.txt';
        $command .= ' -s ' .$sampleName;
        $command .= ' -n ' .LoadConfig::getParam($rH_cfg, 'filterSV', 'minReadSupport');
        $command .= ' -t ' .LoadConfig::getParam($rH_cfg, 'filterSV', 'minReadSupport');
        $command .= ' -b ' .$normalFile
        $command .= ' -c ' .$tumorFile
        $command .= ' && ';
        $command .= filterResults($rH_cfg,  'filterSV', $outputPrefix, '') ;
        $command .= ' && ';
        $command .= generateBedResults($rH_cfg, $outputPrefix) ;
    #}
    return $command;
}


sub filterPI {
    my $rH_cfg          = shift;
    my $sampleName      = shift;
    my $inputPICalls  = shift;
    my $outputPrefix    = shift;
    my $normalFile    = shift;
    my $tumorFile    = shift;

    my $outDate = -M $outputPrefix.'.filteredSV.txt';
    my $inDate = -M $inputDNACCalls;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'filterSV', 'moduleVersion.svtools').' ;';
        $command .= ' \${SVTOOLS_HOME}/Cancer/filterOutPI.py';
        $command .= ' -f ' .$inputPICalls;
        $command .= ' -o ' .$outputPrefix.'.txt';
        $command .= ' -s ' .$sampleName;
        $command .= ' -n ' .LoadConfig::getParam($rH_cfg, 'filterSV', 'minReadSupport');
        $command .= ' -t ' .LoadConfig::getParam($rH_cfg, 'filterSV', 'minReadSupport');
        $command .= ' && ';
        $command .= filterResults($rH_cfg,  'filterSV', $outputPrefix, '') ;
        $command .= ' && ';
        $command .= generateBedResults($rH_cfg, $outputPrefix) ;
    #}
    return $command;
}



sub filterResults {
    my $rH_cfg          = shift;
    my $stepIniPrefix  = shift;
    my $outputPrefix    = shift;
    my $cnvProx           = shift;

    my $outDate = -M $outputPrefix.'.filteredSV.txt';
    my $inDate = -M $inputDNACCalls;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'moduleVersion.svtools').' ;';
        $command .= ' \${SVTOOLS_HOME}/Cancer/filterBedResults.sh';
        $command .= ' '.$outputPrefix.'.txt';
        $command .= ' '.LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceMappabilityBed');
        $command .= ' '.LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceGeneCoordinates');
        $command .= ' '.LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceDGVCoordinates');
        $command .= ' '.LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceMicrosatellitesCoordinates');
        $command .= ' '.LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceRepeatMaskerCoordinates');
        $command .= ' '.LoadConfig::getParam($rH_cfg, $stepIniPrefix, 'referenceGenomeLengths');
        $command .= ' '.$outputPrefix.'.bed';
        $command .= ' '.$outputPrefix.'.tmp';
        $command .= ' '.$cnvProx;
    #}
    return $command;
}

sub generateBedResults {
    my $rH_cfg          = shift;
    my $outputPrefix    = shift;


    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'bedSV', 'moduleVersion.svtools').' ;';
        $command .= ' \${SVTOOLS_HOME}/Cancer/rtxt2rbed.sh ';
        $command .= ' '.$outputPrefix.'.bed.other.filteredSV.annotate.txt';
        $command .= ' '.$outputPrefix.'.bed.other.filteredSV.annotate.bed';
        $command .= ' && '.$outputPrefix.'.bed.TumS.filteredSV.annotate.txt';
        $command .= ' '.$outputPrefix.'.bed.TumS.filteredSV.annotate.bed';
    #}
    return $command;
}


1;
