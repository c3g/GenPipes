#!/usr/env/perl

=head1 NAME

I<MACS2>

=head1 SYNOPSIS

MACS2->getDesign()
MACS2->generateMacs()

=head1 DESCRIPTION

B<MACS2> is a software to execute Model-based Analysis of ChIP-Seq (MACS) on short reads sequencers such as Genome Analyzer (Illumina / Solexa)

=head1 AUTHOR

B<Johanna Sandoval> - I<johanna.sandoval@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package MACS2;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------

sub getDesign {
  my $rH_cfg         = shift;
  my $designFilePath = shift;
  my $designName;
  my $designType;
  my @designComposition;
  my %design;

  open(INFO, $designFilePath) or die "Can't find file $designFilePath\n";
  my @infos = <INFO>;
  close(INFO);
  chomp($infos[0]);
  my @splitA = split(/\t/, $infos[0]);
  my $numberDesigns = scalar(@splitA) - 1;
  for (my $i = 1; $i <= $numberDesigns; $i++) {
    @designComposition = split(",", $splitA[$i]);
    if (scalar(@designComposition) > 1) {
      $designName = $designComposition[0];
      $designType = $designComposition[1];
      if ($designType ne 'B' && $designType ne 'N') {
        die "ERROR: Wrong peak type, should be B (broad) or N (narrow)\n";
      }
    } else {
      die "ERROR: Wrong group definition (name, type); check design file\n";
    }
    my @group1;
    my @group2;
    my @type;
    push(@type, $designType);
    for (my $j = 1; $j < @infos; $j++) {
      chomp($infos[$j]);
      my @splitB = split(/\t/, $infos[$j]);
      my $sampleName = $splitB[0];
      #chomp($sampleName);
      if ($splitB[$i] eq "") {
        die "ERROR: Wrong assignment for a sample; check design file\n";
      } elsif ($splitB[$i] == 1) {
        push(@group1, $sampleName);
      } elsif ($splitB[$i] == 2) {
        push(@group2, $sampleName);
      } elsif ($splitB[$i] == 0) {
        ;  # do nothing
      } else {
        die "ERROR: Wrong group assignment; check design file\n";
      }
    }
    $design{$designName} = [\@group1, \@group2, \@type];
  }
  return \%design;
}


sub generatePeaks {
  my $rH_cfg     = shift;
  my $designName = shift;
  my $treatment  = shift;
  my $control    = shift;
  my $type       = shift;
  my $outputDir  = shift;
  my $paired     = shift;
  my $command    = "";

  my @inputs;
  push(@inputs, $treatment);
  # Compute Genome size or retrieve from config
  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'genomeName');
  my $genomeSize = '';
  if ($refGenome eq 'hg1k_v37' || $refGenome eq 'hg19' || $refGenome eq 'hg18') {
    $genomeSize = ' -g hs ';
  }
  elsif ($refGenome eq 'mm9' || $refGenome eq 'mm10') {
    $genomeSize = ' -g mm ';
  } elsif (LoadConfig::getParam($rH_cfg, 'macs', 'genomeSize') != "") {
    if ($genomeSize > 0){
      $genomeSize = ' -g '.LoadConfig::getParam($rH_cfg, 'macs', 'genomeSize');
    } else {
      die "ERROR: Undefined genome name $refGenome or wrong definition for genome size variable in configuration file (genomeSize) \n";
    }
  } else {
    die "ERROR: Undefined genome name $refGenome or undefined genome size variable in configuration file (genomeSize) \n";
  }
  # Call peaks command
  $command .= ' module load ' . LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.python') . ' ' . LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.macs') . ' && ';
  my $feoptions = "";
  my $extraFlags = "";

  # additional options for MACS must be set in the configuration (ini) file. options --diag femacs femin festep are currently not functional (2013-06-06).
  # FEMIN and FEMAX are the minimum and maximum fold enrichment to consider, and FESTEP is the interval of fold enrichment.
  if (LoadConfig::getParam($rH_cfg, 'macs', 'femin') eq "" || LoadConfig::getParam($rH_cfg, 'macs', 'femax') eq "" || LoadConfig::getParam($rH_cfg, 'macs', 'festep') eq "") {
      $feoptions='';
  } else {
    if (LoadConfig::getParam($rH_cfg, 'macs', 'femin') >= 0 && LoadConfig::getParam($rH_cfg, 'macs', 'femax') >= 0 && LoadConfig::getParam($rH_cfg, 'macs', 'festep') >= 0) {
      $feoptions = ' --fe-min=' . LoadConfig::getParam($rH_cfg, 'macs', 'femin') . ' --fe-max=' . LoadConfig::getParam($rH_cfg, 'macs', 'femax') . ' --fe-step=' . LoadConfig::getParam($rH_cfg, 'macs', 'festep');
    }
  }
  if (!LoadConfig::getParam($rH_cfg, 'macs', 'extraFlags') eq "") {
    $extraFlags = LoadConfig::getParam($rH_cfg, 'macs', 'extraFlags');
  }
  #Peak calling strategies.
  my $options = "";
  my $typeoptions = "";
  # Specific options for treatment / control pairs
  # Paired / single end reads
  if (defined($paired) && $paired == 1) {
    $typeoptions .= " -f BAMPE ";
  } else {
    $typeoptions .= " -f BAM ";
  }

  if (defined($control) && defined($treatment)) {
    push(@inputs, $control);
    if ($type eq 'B') {
      $options = ' --nomodel --broad ';
    } else {
      $options = " ";
    }

    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'macs', 'macsBin') . ' callpeak ' . '-t ' . $treatment . ' -c ' . $control . ' --name=' . $outputDir . '/' . $designName . $genomeSize . $options . ' ' . $typeoptions . ' ' . $feoptions . ' ' . $extraFlags . ' >& ' . $outputDir . '/' . $designName . '.diag.macs.out';
  } elsif (!defined($treatment)) {
    print $treatment;
    die "ERROR: Something wrong with design; treatment should be assigned\n";
  } elsif (!defined($control) && defined($treatment)) {
    if ($type eq 'B') {
      $options = ' --nomodel --nolambda --broad ';
    } else {
      $options = ' --nomodel --nolambda';
    }
    $command .= ' ' . LoadConfig::getParam($rH_cfg, 'macs', 'macsBin') . ' callpeak ' . '-t ' . $treatment . ' --name=' . $outputDir . '/' . $designName . $genomeSize .  $options . ' ' . $typeoptions . ' ' . $feoptions . ' ' . $extraFlags . ' >& ' . $outputDir . '/' . $designName . '.diag.macs.out';
  } elsif (!defined($control) && !defined($treatment)) {
    die "ERROR: Something wrong with design; treatment and control (if available) should be assigned\n";
  }

  my $ro_job = new Job();
  $ro_job->testInputOutputs(\@inputs, [$outputDir . '/' . $designName. '_peaks.xls']);
  if (!$ro_job->isUp2Date()) {
    $ro_job->addCommand($command);
  }

  return $ro_job;
}

1;
