#!/usr/env/perl

=head1 NAME

I<CountIlluminaBarcode>

=head1 SYNOPSIS

CountIlluminaBarcodes->count()

=head1 DESCRIPTION

B<CountIlluminaBarcode> is a library to extract the count of illumina barcodes present in an illumina lane

Output = file

=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package CountIlluminaBarcodes;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub count {
  my $rH_cfg           = shift;
  my $baseCallDir      = shift;
  my $lane             = shift;
  my $mask             = shift;
  my $inputFileToCheck = shift;
  my $outputFile       = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputFileToCheck],[$outputFile]);

  if (!$ro_job->isUp2Date()) {
    my $command;
    
    $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'generateIndexCount','moduleVersion.java') . ' &&'; 
    $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'generateIndexCount', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'generateIndexCount', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'generateIndexCount', 'maxRam').' -jar ' . LoadConfig::getParam($rH_cfg, 'generateIndexCount', 'jar');
    $command .= ' MAX_MISMATCHES=1 NUM_PROCESSORS=12 BARCODE_FILE=' . LoadConfig::getParam($rH_cfg, 'generateIndexCount', 'barcodeFile');
    $command .= ' BASECALLS_DIR='.$baseCallDir.' LANE='.$lane.' READ_STRUCTURE='.$mask.' METRICS_FILE='.$outputFile;
    $command .= ' TMP_DIR='.LoadConfig::getParam($rH_cfg, 'generateIndexCount', 'tmpDir');

    $ro_job->addCommand($command);
  }
  return $ro_job;
}



1;
