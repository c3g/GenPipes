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


  my $latestFile = -M $inputBAM;
  my $output1= -M $rA_outputBAM->[1];
  my $output2= -M $rA_outputBAM->[2];

  my @command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestFile) || !defined($output1) || !defined($output2) || $latestFile > $output1 || $latestFile > $output2) {
    my $Fcommand = 'module load ' .LoadConfig::getParam($rH_cfg, 'wiggle','moduleVersion.samtools') .' ;';
    $Fcommand .= ' samtools view -h -F 256 -f 81 ' . $inputBAM;
    $Fcommand .= ' > ' .$inputBAM .'tmp1.forward.sam &&';
    $Fcommand .= ' samtools view -h -F 256 -f 161 ' . $inputBAM;
    $Fcommand .= ' > ' .$inputBAM .'tmp2.forward.sam &&';
    $Fcommand .= ' cat ' .$inputBAM .'tmp1.forward.sam';
    $Fcommand .= ' ' .$inputBAM .'tmp2.forward.sam';
    $Fcommand .= ' | samtools view -Sb -';
    $Fcommand .= ' > ' .$inputBAM .'tmp1.forward.bam &&';
    $Fcommand .= ' samtools sort ' .$inputBAM .'tmp1.forward.bam';
    $Fcommand .= ' ' .$rA_outputBAM->[0] .' && ';
    $Fcommand .= ' samtools index ' .$rA_outputBAM->[0] .' && ';
    $Fcommand .= ' rm ' .$inputBAM .'tmp*.forward.*am';
    push(@command,$Fcommand);

    my $Rcommand = 'module load ' .LoadConfig::getParam($rH_cfg, 'wiggle','moduleVersion.samtools') .' ;';
    $Rcommand .= ' samtools view -h -F 256 -f 97 ' . $inputBAM;
    $Rcommand .= ' > ' .$inputBAM .'tmp1.reverse.sam &&';
    $Rcommand .= ' samtools view -h -F 256 -f 145 ' . $inputBAM;
    $Rcommand .= ' > ' .$inputBAM .'tmp2.reverse.sam &&';
    $Rcommand .= ' cat ' .$inputBAM .'tmp1.reverse.sam';
    $Rcommand .= ' ' .$inputBAM .'tmp2.reverse.sam';
    $Rcommand .= ' | samtools view -Sb -';
    $Rcommand .= ' > ' .$inputBAM .'tmp1.reverse.bam &&';
    $Rcommand .= ' samtools sort ' .$inputBAM .'tmp1.reverse.bam';
    $Rcommand .= ' ' .$rA_outputBAM->[1] .' &&';
    $Rcommand .= ' samtools index ' .$rA_outputBAM->[1] .' &&';
    $Rcommand .= ' rm ' .$inputBAM .'tmp*.reverse.*am';
    push(@command,$Rcommand);
  }
    
  return \@command;
}

sub graph{
  my $rH_cfg         = shift;
  my $sampleName     = shift;
  my $inputBAM       = shift;
  my $outputBegGraph = shift;
  my $outputWiggle   = shift;


  my $latestFile = -M $inputBAM;
  my $output1= -M $outputBegGraph;
  my $output2= -M $outputWiggle;

  my $command;
  # -M gives modified date relative to now. The bigger the older.
  if(!defined($latestFile) || !defined($output1) || !defined($output2) || $latestFile > $output1 || $latestFile > $output2) {
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'wiggle','moduleVersion.samtools') .' ' .LoadConfig::getParam($rH_cfg, 'wiggle','moduleVersion.bedtools')  .' ' .LoadConfig::getParam($rH_cfg, 'wiggle','moduleVersion.bed2wig') .' ;';
    $command .= ' nmblines=\$(samtools view -F 256 -f 81 ' . $inputBAM .' | wc -l) &&';
    $command .= ' scalefactor=0\$(echo \"scale=2; 1 / (\$nmblines / 10000000);\" | bc) &&';   
    $command .= ' genomeCoverageBed -bg -ibam ' . $inputBAM;
    $command .= ' -g ' .LoadConfig::getParam($rH_cfg, 'wiggle','chromosomeSizeFile'); 
    $command .= ' -split -scale \$scalefactor ';
    $command .= ' > ' .$outputBegGraph . ' &&';
    $command .= ' bedGraphToBigWig ' .$outputBegGraph;
    $command .= '  ' .LoadConfig::getParam($rH_cfg, 'wiggle','chromosomeSizeFile');
    $command .= '  ' .$outputWiggle;
  }
    
  return $command;
}


 

1;
