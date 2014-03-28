#!/usr/bin/env perl

=head1 NAME

I<ReadSet>

=head1 SYNOPSIS

Object used to hold information on a ReadSet

=head1 DESCRIPTION

Object used to hold information on a ReadSet


=head1 DEPENDENCY

=cut

package ReadSet;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin";

# Dependencies
#-----------------------

# SUB
#-----------------------
sub new {
  my $class = shift;
  my $name = shift;
  my $runType = shift;

  unless ($name =~ /^\w[\w.-]*$/) {
    die "Error: readSet name \"$name\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!";
  }

  unless ($runType eq "PAIRED_END" or $runType eq "SINGLE_END") {
    die "Error: readSet runType \"$runType\" is invalid (should be \"PAIRED_END\" or \"SINGLE_END\")!";
  }

  my $self = {
    '_name' => $name,
    '_runType' => $runType,
  };

  bless($self, $class);
  return $self;
}

sub getName {
  my ($self) = @_;
  return $self->{'_name'};
}

sub getRunType {
  my ($self) = @_;
  return $self->{'_runType'};
}

sub setSample {
  my ($self, $sample) = @_;
  $self->{'_sample'} = $sample;
}

sub getSample {
  my ($self) = @_;
  return $self->{'_sample'};
}

sub setBAM {
  my ($self, $BAM) = @_;
  $self->{'_BAM'} = $BAM;
}

sub getBAM {
  my ($self) = @_;
  return $self->{'_BAM'};
}

sub setFASTQ1 {
  my ($self, $FASTQ1) = @_;
  $self->{'_FASTQ1'} = $FASTQ1;
}

sub getFASTQ1 {
  my ($self) = @_;
  return $self->{'_FASTQ1'};
}

sub setFASTQ2 {
  my ($self, $FASTQ2) = @_;
  $self->{'_FASTQ2'} = $FASTQ2;
}

sub getFASTQ2 {
  my ($self) = @_;
  return $self->{'_FASTQ2'};
}

sub setLibrary {
  my ($self, $library) = @_;
  $self->{'_library'} = $library;
}

sub getLibrary {
  my ($self) = @_;
  return $self->{'_library'};
}

sub setRun {
  my ($self, $run) = @_;
  $self->{'_run'} = $run;
}

sub getRun {
  my ($self) = @_;
  return $self->{'_run'};
}

sub setLane {
  my ($self, $lane) = @_;
  $self->{'_lane'} = $lane;
}

sub getLane {
  my ($self) = @_;
  return $self->{'_lane'};
}

sub setAdaptor1 {
  my ($self, $adaptor1) = @_;
  $self->{'_adaptor1'} = $adaptor1;
}

sub getAdaptor1 {
  my ($self) = @_;
  return $self->{'_adaptor1'};
}

sub setAdaptor2 {
  my ($self, $adaptor2) = @_;
  $self->{'_adaptor2'} = $adaptor2;
}

sub getAdaptor2 {
  my ($self) = @_;
  return $self->{'_adaptor2'};
}

sub setQualityOffset {
  my ($self, $qualityOffset) = @_;
  $self->{'_qualityOffset'} = $qualityOffset;
}

sub getQualityOffset {
  my ($self) = @_;
  return $self->{'_qualityOffset'};
}

sub setBED {
  my ($self, $BED) = @_;
  $self->{'_BED'} = $BED;
}

sub getBED {
  my ($self) = @_;
  return $self->{'_BED'};
}

1;
