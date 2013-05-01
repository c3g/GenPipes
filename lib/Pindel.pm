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
use File::Basename;

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
    my $inDate = -M $tumorBam;
    my $inDate2 = -M $normalMetrics;
    my $inDate3 = -M $tumorMetrics;
    my $inDate4 = -M $normalBam;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    if(!defined($outDate) || !defined($inDate) || $inDate < $outDate || !defined($inDate2) || $inDate2 < $outDate || !defined($inDate3) || $inDate3 < $outDate || !defined($inDate4) || $inDate4 < $outDate) {
        $command .= 'TinS=\$(grep -A 1 \"MEDIAN\"  ' .$tumorMetrics .' | awk \' NR == 2 {print \$1} \') && ';
        $command .= 'NinS=\$(grep -A 1 \"MEDIAN\"  ' .$normalMetrics .' | awk \' NR == 2 {print \$1} \') && ';
        $command .= 'echo -e \"' .$tumorBam .'\t\${TinS}\tTUMOR\n' .$normalBam .'\t\${NinS}\tBLOOD\n\" > ' . $output ;
    }
    
    my $bai1Date = -M $tumorBam .'.bai';
    if(!defined($bai1Date) ) {
        (my $tumorBai = $tumorBam) =~ s/\.bam/\.bai/g ;
        if(defined($command) && length($command) > 0) {
            $command .= ' &&';
        }
        $command .= ' ln -s ' .basename($tumorBai) .' ' .$tumorBam .'.bai' ;
    }
    my $bai2Date = -M $normalBam .'.bai';
    if(!defined($bai2Date) ) {
        if(defined($command) && length($command) > 0) {
            $command .= ' &&';
        }
        (my $normalBai = $normalBam) =~ s/\.bam/\.bai/g ;
        $command .= ' ln -s ' .basename($normalBai) .' ' .$normalBam .'.bai' ;
    }
    return $command;
}


sub pairedPI {
    my $rH_cfg          = shift;
    my $chr             = shift;
    my $inputCFG        = shift;
    my $outputPrefix    = shift;
    my $outputTest     = shift;
    my $PIOption       = shift;

    my $outDate = -M $outputTest .'_SI';
    my $inDate = -M $inputCFG;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
      $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pindel', 'moduleVersion.pindel').' ;';
      $command .= ' pindel '.LoadConfig::getParam($rH_cfg, 'pindel', 'piParameters');
      $command .= ' -f '.$chr;
      $command .= ' -i '.$inputCFG;
      $command .= ' -o '.$outputPrefix;
      $command .= ' '.$PIOption;

    }
    return $command;
}

sub mergeChro {
    my $rH_cfg          = shift;
    my $outputPrefix    = shift;

    my $outDate = -M $outputPrefix .'_SI';
    my $inDate = -M $outputPrefix .'.1_SI';

    my $command;
    # -M gives modified date relative to now. The bigger the older.
    if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
      $command .= 'rm ' .$outputPrefix .'_BP' .' ' .$outputPrefix .'_D' .' ' .$outputPrefix .'_INV' .' ' .$outputPrefix .'_LI' .' ' .$outputPrefix .'_SI' .' ' .$outputPrefix .'_TD' .' ; ' ;
      $command .= 'touch ' .$outputPrefix .'_BP' .' ' .$outputPrefix .'_D' .' ' .$outputPrefix .'_INV' .' ' .$outputPrefix .'_LI' .' ' .$outputPrefix .'_SI' .' ' .$outputPrefix .'_TD' .' && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_BP ; do cat \$i >> '  .$outputPrefix .'_BP  ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_D ; do cat \$i >> '  .$outputPrefix .'_D ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_INV ; do cat \$i >> '  .$outputPrefix .'_INV ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_LI ; do cat \$i >> '  .$outputPrefix .'_LI ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_SI ; do cat \$i >> '  .$outputPrefix .'_SI ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_TD ; do cat \$i >> '  .$outputPrefix .'_TD ; done ' ;
    }
    return $command;
}

    
1;
