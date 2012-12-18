#!/usr/env/perl

=head1 NAME

I<IGVTools>

=head1 SYNOPSIS

IGVTools::computeTDF

=head1 DESCRIPTION

B<IGVTools> is a library used to generate tiled files for IGV

Input = file_name

Output = array

=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package IGVTools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub computeTDF {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $inputBAM    = shift;
# DOES NOT WORK YET!!!

  my $command;
  $command .= 'module load mugqic/igvtools/2.1.24 ;';
  $command .= ' igvtools count -f min,max,mean ';
  $command .= $inputBAM.' '.$inputBAM.'.tdf';
  $command .= ' '.LoadConfig::getParam($rH_cfg, 'computeTDF', 'igvGenome');
  return $command;
}

1;
