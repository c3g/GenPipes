#!/usr/env/perl

=head1 NAME

I<Wiggle>

=head1 SYNOPSIS

Wiggle-> strandBam()
Wiggle-> graph()

=head1 DESCRIPTION

B<Wiggle> is a library to generate wiggle graphs

Input = file_name

Output = array


=head1 AUTHOR
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package Wiggle;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------
sub strandBam{
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $inputBAM      = shift;
  my $rA_outputBAM  = shift;


  my $latestFile = -M $inputFile;
  my $output1= -M $rA_outputBAM[1];
  my $output2= -M $rA_outputBAM[2];

  my @command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestFile) || !defined($output1) || !defined($output2) || $latestFile > $output1 || $latestFile > $output2) {
    my $Fcommand = 'module load ' .LoadConfig::getParam($rH_cfg, 'wiggle','samtoolsModule') .' ;';
    $Fcommand .= ' samtools view -h -F 256 -f 81 ' . $inputBAM;
    $Fcommand .= ' > alignment/' .$ampleName .'/' .$ampleName .'.tmp1.forward.sam ;';
    $Fcommand .= ' samtools view -h -F 256 -f 161 ' . $inputBAM;
    $Fcommand .= ' > alignment/' .$ampleName .'/' .$ampleName .'.tmp2.forward.sam ;';
    $Fcommand .= ' cat alignment/' .$ampleName .'/' .$ampleName .'.tmp1.forward.sam';
    $Fcommand .= ' alignment/' .$ampleName .'/' .$ampleName .'.tmp2.forward.sam';
    $Fcommand .= ' | samtools view -Sb -';
    $Fcommand .= ' > alignment/' .$ampleName .'/' .$ampleName .'.tmp1.forward.bam ;';
    $Fcommand .= ' samtools sort alignment/' .$ampleName .'/' .$ampleName .'.tmp1.forward.bam';
    $Fcommand .= ' ' .$rA_outputBAM[0] .' ; ';
    $Fcommand .= ' samtools index ' .$rA_outputBAM->[0] .' ; ';
    $Fcommand .= ' rm alignment/' .$ampleName .'/' .$ampleName .'.tmp*.forward.*am';
    push(@command,$Fcommand);

    my $Rcommand = 'module load ' .LoadConfig::getParam($rH_cfg, 'wiggle','samtoolsModule') .' ;';
    $Rcommand .= ' samtools view -h -F 256 -f 97 ' . $inputBAM;
    $Rcommand .= ' > alignment/' .$ampleName .'/' .$ampleName .'.tmp1.reverse.sam ;';
    $Rcommand .= ' samtools view -h -F 256 -f 145 ' . $inputBAM;
    $Rcommand .= ' > alignment/' .$ampleName .'/' .$ampleName .'.tmp2.reverse.sam ;';
    $Rcommand .= ' cat alignment/' .$ampleName .'/' .$ampleName .'.tmp1.reverse.sam';
    $Rcommand .= ' alignment/' .$ampleName .'/' .$ampleName .'.tmp2.reverse.sam';
    $Rcommand .= ' | samtools view -Sb -';
    $Rcommand .= ' > alignment/' .$ampleName .'/' .$ampleName .'.tmp1.reverse.bam ;';
    $Rcommand .= ' samtools sort alignment/' .$ampleName .'/' .$ampleName .'.tmp1.reverse.bam';
    $Rcommand .= ' ' .$rA_outputBAM[1] .' ; ';
    $Rcommand .= ' samtools index ' .$rA_outputBAM->[1] .' ; ';
    $Rcommand .= ' rm alignment/' .$ampleName .'/' .$ampleName .'.tmp*.reverse.*am';
    push(@command,$Rcommand);
  }
    
  return \@commands;
}

sub graph{
  my $rH_cfg         = shift;
  my $sampleName     = shift;
  my $inputBAM       = shift;
  my $outputBegGraph = shift;
  my $outputWiggle   = shift;


  my $latestFile = -M $inputFile;
  my $output1= -M $outputBegGraph;
  my $output2= -M $outputWiggle;

  my $command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestFile) || !defined($output1) || !defined($output2) || $latestFile > $output1 || $latestFile > $output2) {
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'wiggle','samtoolsModule') .' ' .LoadConfig::getParam($rH_cfg, 'wiggle','bedtoolsModule')  .' ' .LoadConfig::getParam($rH_cfg, 'wiggle','bed2wigModule') .' ;';
    $command .= ' nmblines=\$(samtools view -F 256 -f 81 ' . $inputBAM .' | wc -l) ;';
    $command .= ' scalefactor=0\$(echo \"scale=2; 1 / (\$nmblines / 10000000);\" | bc) ;';   
    $command .= ' genomeCoverageBed -bg -ibam ' . $inputBAM;
    $command .= ' -g ' .LoadConfig::getParam($rH_cfg, 'wiggle','chromosomeSizeFile'); 
    $command .= ' -split -scale \$scalefactor ';
    $command .= ' > ' .$outputBegGraph . ' ;';
    $command .= ' bedGraphToBigWig ' .$outputBegGraph;
    $command .= '  ' .LoadConfig::getParam($rH_cfg, 'wiggle','chromosomeSizeFile');
    $command .= '  ' .$outputWiggle;
  }
    
  return $commands;
}


 

1;
