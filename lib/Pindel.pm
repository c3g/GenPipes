#!/usr/env/perl

=head1 NAME

I<Pindel>

=head1 SYNOPSIS

Pindel

=head1 DESCRIPTION

B<Pindelr> is a library to analyse SV events in a genome

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package Pindel;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub pairedConfigFile {
    my $rH_cfg        = shift;
    my $tumorMetrics  = shift;
    my $normalMetrics = shift;
    my $tumorBam      = shift;
    my $normalBam      = shift;
    my $output        = shift;

    my $outDate = -M $output;
    my $inDate = -M $tumorBAM;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'TinS=\$(grep -A 1 \"MEDIAN\"  ' .$tumorMetrics .' | awk \' NR == 2 {print \$1} \')\' && ';
        $command .= 'NinS=\$(grep -A 1 \"MEDIAN\"  ' .$normalMetrics .' | awk \' NR == 2 {print \$1} \')\' && ';
        $command .= 'echo -e \"' .$tumorBam .'\t\${TinS}\tTUMOR\n' .$normalBam .'\t\${NinS}\tBLOOD\n\" > ' . $output ;
    #}
    return $command;
}


sub pairedPI {
    my $rH_cfg          = shift;
    my $chr             = shift;
    my $inputCFG        = shift;
    my $outputPrefix    = shift;
    my $brdOption       = shift;

#     my $outDate = -M $outputPrefix.'.ctx';
#     my $inDate = -M $inputCFG;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
    $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pindel', 'moduleVersion.pindel').' ;';
    $command .= ' pindel '.LoadConfig::getParam($rH_cfg, 'pindel', 'piParameters');
    $command .= ' -f '.$chr;
    $command .= ' -i '.$inputCFG;
    $command .= ' -o '.$outputPrefix;
    $command .= ' '.$brdOption;

    #}
    return $command;
}

sub mergeChro {
    my $rH_cfg          = shift;
    my $outputPrefix    = shift;

    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
    $command .= 'rm ' .$outputPrefix .'_BP' .' ' .$outputPrefix .'_D' .' ' .$outputPrefix .'_INV' .' ' .$outputPrefix .'_LI' .' ' .$outputPrefix .'_SI' .' ' .$outputPrefix .'_TD' .' && ' ;
    $command .= 'touch ' .$outputPrefix .'_BP' .' ' .$outputPrefix .'_D' .' ' .$outputPrefix .'_INV' .' ' .$outputPrefix .'_LI' .' ' .$outputPrefix .'_SI' .' ' .$outputPrefix .'_TD' .' && ' ;
    $command .= 'for i in ' .$outputPrefix .'.*_BP ; do cat \$i >> '  .$outputPrefix .'_BP && ' ;
    $command .= 'for i in ' .$outputPrefix .'.*_D ; do cat \$i >> '  .$outputPrefix .'_D && ' ;
    $command .= 'for i in ' .$outputPrefix .'.*_INV ; do cat \$i >> '  .$outputPrefix .'_INV && ' ;
    $command .= 'for i in ' .$outputPrefix .'.*_LI ; do cat \$i >> '  .$outputPrefix .'_LI && ' ;
    $command .= 'for i in ' .$outputPrefix .'.*_SI ; do cat \$i >> '  .$outputPrefix .'_SI && ' ;
    $command .= 'for i in ' .$outputPrefix .'.*_TD ; do cat \$i >> '  .$outputPrefix .'_TD' ;
    #}
    return $command;
}

    
1;
