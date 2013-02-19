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
  
    my $outDate = -M $output;
    my $inDate = -M $sampleBAM;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'bam2cfg', 'moduleVersion.breakdancer').' ;';
        $command .= ' module load '.LoadConfig::getParam($rH_cfg, 'bam2cfg', 'moduleVersion.samtools').' ;';
        $command .= ' bam2cfg.pl -g -h ';
        $command .= ' -c '.$stdDevCutoff;
        $command .= ' '.$sampleBAM;
        $command .= ' > '.$output;
    #}
    return $command;
}

sub pairedBRDITX {
    my $rH_cfg          = shift;
    my $sampleName      = shift;
    my $inputCFG        = shift;
    my $outputPrefix    = shift;

    my $outDate = -M $outputPrefix.'.ctx';
    my $inDate = -M $inputCFG;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pairedBRDITX', 'moduleVersion.breakdancer').' ;';
        $command .= ' breakdancer_max '.LoadConfig::getParam($rH_cfg, 'pairedBRDITX', 'brdParameters');
        $command .= ' -g '.$outputPrefix.'.bed';
        $command .= ' -d '.$outputPrefix.'.ctx';
        $command .= ' '.$inputCFG;
        $command .= ' > '.$outputPrefix.'.ctx';
    #}
    return $command;
}

sub pairedBRD {
    my $rH_cfg          = shift;
    my $sampleName      = shift;
    my $chr             = shift;
    my $inputCFG        = shift;
    my $outputPrefix    = shift;

    my $outDate = -M $outputPrefix.'.ctx';
    my $inDate = -M $inputCFG;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
    $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pairedBRD', 'moduleVersion.breakdancer').' ;';
    $command .= ' breakdancer_max '.LoadConfig::getParam($rH_cfg, 'pairedBRD', 'brdParameters');
    $command .= ' -o '.$chr;
    $command .= ' -g '.$outputPrefix.'.bed';
    $command .= ' -d '.$outputPrefix.'.ctx';
    $command .= ' '.$inputCFG;
    $command .= ' > '.$outputPrefix.'.ctx';
    #}
    return $command;
}
1;
