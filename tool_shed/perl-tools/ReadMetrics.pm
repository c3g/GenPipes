#!/usr/env/perl

=head1 NAME

I<ReadMetrics>

=head1 SYNOPSIS

ReadMetrics-> parseFlagstats()
ReadMetrics-> parseTrimOutput()
ReadMetrics-> mergeStats()

=head1 DESCRIPTION

B<ReadMetrics> is a library to generate QC, stats and metrics


=head1 AUTHOR
B<Johanna Sandoval> - I<johanna.sandoval@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package ReadMetrics;

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

sub parseFlagstats{
	my $sampleName    = shift;
	my $flagStatsFile = shift;
  my %stats         ;
  my $nbQCPassedReads = -1;
  my $nbDuplicateReads= -1;
  my $nbAlignedReads  = -1;
  open(FLAGSTATS, "$flagStatsFile") or die "Can't open $flagStatsFile\n";
  while(my $line = <FLAGSTATS>) {
    chomp($line);
    if($line =~ /^([0-9]+) \+ [0-9]+ in total/) {
       $nbQCPassedReads = $1+0;
    }
    elsif($line =~ /^([0-9]+) \+ [0-9]+ duplicates/) {
      $nbDuplicateReads += $1;
    }
    elsif($line =~ /^([0-9]+) \+ [0-9]+ mapped/) {
      $nbAlignedReads += $1;
    }
  }
  $stats{$sampleName}=[$nbQCPassedReads, $nbAlignedReads, $nbDuplicateReads];
  close(FLAGSTATS);
  return \%stats;
}

sub parseTrimOutput{
  my $sampleName      = shift;
  my $trimOutputFile  = shift;
  my %stats           ;
  my $nbRawReads      = -1;
  my $nbFilteredReads = -1;
  my $nbSingleFiltered= -1;
  
  open(TRIMSTATS, "$trimOutputFile") or die "Can't open $trimOutputFile\n";
  while(my $line = <TRIMSTATS>) {
    chomp($line);
    if($line =~ /^Raw Fragments\,([0-9]+)/){
      $nbRawReads = $1;
    }elsif($line =~ /^Fragment Surviving\,([0-9]+)/ ){
      $nbFilteredReads = $1;
    }elsif($line =~ /^Single Surviving\,([0-9]+)/ ){
      $nbSingleFiltered= $1;
    }    
  }
  $stats{$sampleName}=[$nbRawReads, $nbFilteredReads, $nbSingleFiltered];
  close(TRIMSTATS);  
  return \%stats;
}

sub mergeStats{
  my $sampleName           = shift;
  my $outputFile           = shift;
  my $rHoA_trimStats       = shift;
  my $rHoA_alignmentStats  = shift;
  
  # Open an output file
  open(OUTFILE, ">".$outputFile);
	print OUTFILE $sampleName.",".join(',', @{$rHoA_trimStats->{$sampleName}}). ','. join(',', @{$rHoA_alignmentStats->{$sampleName}})."\n";
	close(OUTFILE);
}

1;