#!/usr/env/perl

=head1 NAME

I<BWA>

=head1 SYNOPSIS

BWA->aln()

=head1 DESCRIPTION

B<BWA> is a library that aligns fastqs on a reference genome

Input = file_name

Output = array


=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package BWA;

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
sub aln {
  my $rH_cfg        = shift;
  my $sampleName    = shift;
  my $rH_laneInfo   = shift;
  my $pair1         = shift;
  my $pair2         = shift;
  my $single1       = shift;
  my $single2       = shift;
  my $optOutputTag  = shift;

  my $command = "";

  if ( $rH_laneInfo->{'runType'} eq "SINGLE_END" ) {
    $command = singleCommand($rH_cfg, $sampleName, $rH_laneInfo, $single1, $optOutputTag);
  }
  elsif($rH_laneInfo->{'runType'} eq "PAIRED_END") {
    $command = pairCommand($rH_cfg, $sampleName, $rH_laneInfo, $pair1, $pair2, $optOutputTag);
  }
  else {
    die "Unknown runType: ".$rH_laneInfo->{' runType '}."\n";
  }
    
  return $command;
}

sub pairCommand {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rH_laneInfo = shift;
  my $pair1       = shift;
  my $pair2       = shift;
  my $optOutputTag= shift;

  my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
  my $outputSai1Name = $laneDirectory . $sampleName.'.pair1.sai';
  my $outputSai2Name = $laneDirectory . $sampleName.'.pair2.sai';
  my $outputBAM = $laneDirectory . $sampleName.$optOutputTag.'.sorted.bam';
  my $bamFileDate = -M $outputBAM;

  my @commands;
  # -M gives modified date relative to now. The bigger the older.
  if (!defined($bamFileDate) || !defined(-M $pair1) || !defined(-M $pair2) || $bamFileDate > -M $pair1 || $bamFileDate > -M $pair2) {
    my $sai1Command = "";
    my $sai2Command = "";
    my $bwaCommand = "";

    $sai1Command .= 'module load mugqic/bwa/0.6.2 ; bwa aln';
    $sai1Command .= ' -t '.LoadConfig::getParam($rH_cfg, 'aln', 'bwaAlnThreads');
    $sai1Command .= ' '.LoadConfig::getParam($rH_cfg, 'aln', 'bwaRefIndex');
    $sai1Command .= ' '.$pair1;
    $sai1Command .= ' -f '.$outputSai1Name;
    push(@commands, $sai1Command);
    
    $sai2Command .= 'module load mugqic/bwa/0.6.2 ; bwa aln';
    $sai2Command .= ' -t '.LoadConfig::getParam($rH_cfg, 'aln', 'bwaAlnThreads');
    $sai2Command .= ' '.LoadConfig::getParam($rH_cfg, 'aln', 'bwaRefIndex');
    $sai2Command .= ' '.$pair2;
    $sai2Command .= ' -f '.$outputSai2Name;
    push(@commands, $sai2Command);
    
    my $rgId = $rH_laneInfo->{'libraryBarcode'}."_".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'};
    my $rgTag = '\"@RG\tID:'.$rgId.'\tSM:'.$rH_laneInfo->{'name'}.'\tLB:'.$rH_laneInfo->{'libraryBarcode'}.'\tPU:run'.$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'}.'\tCN:'.LoadConfig::getParam($rH_cfg, 'aln', 'bwaInstitution').'\tPL:Illumina\"';
    $bwaCommand .= 'module load mugqic/bwa/0.6.2 ;module load mugqic/picard/1.77 ; bwa sampe';
    $bwaCommand .= ' -r '.$rgTag;
    $bwaCommand .= ' '.LoadConfig::getParam($rH_cfg, 'aln', 'bwaRefIndex');
    $bwaCommand .= ' '.$outputSai1Name;
    $bwaCommand .= ' '.$outputSai2Name;
    $bwaCommand .= ' '.$pair1;
    $bwaCommand .= ' '.$pair2;
    $bwaCommand .= ' | java -Xmx'.LoadConfig::getParam($rH_cfg, 'aln', 'sortSamRam').' -jar \${PICARD_HOME}/SortSam.jar INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate';
    $bwaCommand .= ' OUTPUT='.$outputBAM;
    $bwaCommand .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'aln', 'sortSamRecInRam');
    push(@commands, $bwaCommand);
  }

  return \@commands;
}

sub singleCommand {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $rH_laneInfo = shift;
  my $single      = shift;
  my $optOutputTag= shift;

  my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
  my $outputSaiName = $laneDirectory . $sampleName.'.single.sai';
  my $outputBAM = $laneDirectory . $sampleName.$optOutputTag.'.sorted.bam';
  my $bamFileDate = -M $outputBAM;

  my @commands;
  # -M gives modified date relative to now. The bigger the older.
  if ($bamFileDate > -M $single) {
    my $saiCommand = "";
    my $bwaCommand = "";

    $saiCommand .= 'module load mugqic/bwa/0.6.2 ; bwa aln';
    $saiCommand .= ' -t '.LoadConfig::getParam($rH_cfg, 'aln', 'bwaAlnThreads');
    $saiCommand .= ' '.LoadConfig::getParam($rH_cfg, 'aln', 'bwaRefIndex');
    $saiCommand .= ' '.$single;
    $saiCommand .= ' -f '.$outputSaiName;
    push(@commands, $saiCommand);
    
    my $rgId = $rH_laneInfo->{'libraryBarcode'}."_".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'};
    my $rgTag = '\"@RG\tID:'.$rgId.'\tSM:'.$rH_laneInfo->{'name'}.'\tLB:'.$rH_laneInfo->{'libraryBarcode'}.'\tPU:run'.$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'}.'\tCN:'.LoadConfig::getParam($rH_cfg, 'aln', 'bwaInstitution').'\tPL:Illumina\"';
    $bwaCommand .= 'module load mugqic/bwa/0.6.2 ;module load mugqic/picard/1.77 ; bwa samse';
    $bwaCommand .= ' -r '.$rgTag;
    $bwaCommand .= ' '.LoadConfig::getParam($rH_cfg, 'aln', 'bwaRefIndex');
    $bwaCommand .= ' '.$outputSaiName;
    $bwaCommand .= ' '.$single;
    $bwaCommand .= ' | java -Xmx'.LoadConfig::getParam($rH_cfg, 'aln', 'sortSamRam').' -jar \${PICARD_HOME}/SortSam.jar INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate';
    $bwaCommand .= ' OUTPUT='.$outputBAM;
    $bwaCommand .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'aln', 'sortSamRecInRam');
    push(@commands, $bwaCommand);
  }

  return \@commands;
}

sub index{
	my $rH_cfg      = shift;
	my $sampleName  = shift;
	my $rH_laneInfo = shift;
	
	my $laneDirectory = $sampleName . "/run" . $rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} . "/";
    my %retVal;
    
    my $command = 'module add mugqic/bwa/0.6.2;';
    $command .= ' bwa index ' . $laneDirectory . $rH_cfg->{'aln.bwaRefIndex'};
    
    $retVal{'command'} = $command;
    return (\%retVal);
    
}



1;











