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
use PipelineUtils;
use File::Basename;

# SUB
#-----------------------
sub pairedConfigFile {
    my $rH_cfg        = shift;
    my $tumorMetrics  = shift;
    my $normalMetrics = shift;
    my $tumorBam      = shift;
    my $normalBam      = shift;
    my $output        = shift;

    my $up2date = PipelineUtils::testInputOutputs([$tumorBam,$normalMetrics,$tumorMetrics,$normalBam], [$output]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'TinS=\$(grep -A 1 \"MEDIAN\"  ' .$tumorMetrics .' | awk \' NR == 2 {print \$1} \') && ';
        $command .= 'NinS=\$(grep -A 1 \"MEDIAN\"  ' .$normalMetrics .' | awk \' NR == 2 {print \$1} \') && ';
        $command .= 'echo -e \"' .$tumorBam .'\t\${TinS}\tTUMOR\n' .$normalBam .'\t\${NinS}\tBLOOD\n\" > ' . $output ;
    
        my ($tumorBai) = $tumorBam =~ s/\.bam/\.bai/g ;
        $command .= ' && rm -f '.$tumorBam .'.bai && ln -s ' .basename($tumorBai) .' ' .$tumorBam .'.bai' ;
        my ($normalBai) = $normalBam =~ s/\.bam/\.bai/g ;
        $command .= ' && rm -f '.$normalBam .'.bai && ln -s ' .basename($normalBai) .' ' .$normalBam .'.bai' ;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }
    return $ro_job;
}


sub pairedPI {
    my $rH_cfg          = shift;
    my $chr             = shift;
    my $inputCFG        = shift;
    my $outputPrefix    = shift;
    my $outputTest     = shift;
    my $PIOption       = shift;

    my $up2date = PipelineUtils::testInputOutputs([$inputCFG], [$outputTest .'_SI']);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'pindel', 'moduleVersion.pindel').' ;';
      $command .= ' pindel '.LoadConfig::getParam($rH_cfg, 'pindel', 'piParameters');
      $command .= ' -f '.$chr;
      $command .= ' -i '.$inputCFG;
      $command .= ' -o '.$outputPrefix;
      $command .= ' '.$PIOption;
      $command .= ' ' . $up2date;

      $ro_job->addCommand($command);

    }
    return $ro_job;
}

sub mergeChro {
    my $rH_cfg          = shift;
    my $outputPrefix    = shift;

    my $up2date = PipelineUtils::testInputOutputs([$outputPrefix .'.1_SI'], [$outputPrefix .'_SI']);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
      my $command;
      $command .= 'rm ' .$outputPrefix .'_BP' .' ' .$outputPrefix .'_D' .' ' .$outputPrefix .'_INV' .' ' .$outputPrefix .'_LI' .' ' .$outputPrefix .'_SI' .' ' .$outputPrefix .'_TD' .' ; ' ;
      $command .= 'touch ' .$outputPrefix .'_BP' .' ' .$outputPrefix .'_D' .' ' .$outputPrefix .'_INV' .' ' .$outputPrefix .'_LI' .' ' .$outputPrefix .'_SI' .' ' .$outputPrefix .'_TD' .' && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_BP ; do cat \$i >> '  .$outputPrefix .'_BP  ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_D ; do cat \$i >> '  .$outputPrefix .'_D ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_INV ; do cat \$i >> '  .$outputPrefix .'_INV ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_LI ; do cat \$i >> '  .$outputPrefix .'_LI ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_SI ; do cat \$i >> '  .$outputPrefix .'_SI ; done && ' ;
      $command .= 'for i in ' .$outputPrefix .'.*_TD ; do cat \$i >> '  .$outputPrefix .'_TD ; done ' ;
      $command .= ' ' . $up2date;

      $ro_job->addCommand($command);
    }
    return $ro_job;
}
    
1;
