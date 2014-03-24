#!/usr/bin/env perl

=head1 NAME

I<Pipeline>

=head1 SYNOPSIS

Object used to hold information on a Pipeline to run

=head1 DESCRIPTION

Object used to hold information on a Pipeline to run


=head1 DEPENDENCY

=cut

package Pipeline;

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
use Cwd 'abs_path';
use File::Basename;
use Parse::Range qw(parse_range);
use LoadConfig;
use Sample;
use Step;

# SUB
#-----------------------
sub new {
  my $class = shift;
  my $rA_steps = shift;
  my $sampleFile = shift;
  my $workDir = shift;

  my $self = {
    '_steps' => [],
    '_samples' => [],
  };

  bless($self, $class);

  for my $rH_step (@$rA_steps) {
    # Retrieve parent step objects from names; raise error if they don't exist
    my $parentSteps = [];

    foreach my $parentStepName (@{$rH_step->{'parentSteps'}}) {
      my $parentStep = $self->getStepByName($parentStepName);

      if ($parentStep) {
        push(@$parentSteps, $parentStep);
      } else {
        die "Error: parent step \"$parentStepName\" of step \"" . $rH_step->{'name'} . "\" does not exist!";
      }
    }

    $self->addStep(Step->new($rH_step->{'name'}, $rH_step->{'loop'}, $parentSteps));
  }

  # Create and add samples to this pipeline
  foreach my $sample (@{Sample->parseSampleFile($sampleFile)}) {
    $self->addSample($sample);
  }

  # Check working directory and set it to current one if not defined
  if (defined $workDir) {
    if (-d $workDir) {
      $workDir = abs_path($workDir);
    } else {
      die "Error: working directory $workDir does not exist or is not a directory";
    }
  } else {
    $workDir = "`pwd`";
  }

  # Add script name (without suffix) as job list filename prefix (in practice, identical to pipeline name)
  my $jobListPrefix = fileparse($0, qr/\.[^.]*/);

  # Set pipeline header and global variables
  print <<END;
#!/bin/bash

WORK_DIR=$workDir
JOB_OUTPUT_ROOT=\$WORK_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=\$JOB_OUTPUT_ROOT/${jobListPrefix}_job_list_\$TIMESTAMP
cd \$WORK_DIR

END

  return $self;
}

sub getSteps {
  my ($self) = @_;
  return $self->{'_steps'};
}

sub getStepsByRange {
  my ($self, $range) = @_;
  # Add undef element to shift 1st step at position 1 instead of 0
  return (undef, @{$self->getSteps()})[parse_range($range)];
}

sub getNbSteps {
  my ($self) = @_;
  return scalar(@{$self->{'_steps'}});
}

sub addStep {
  my ($self, $rO_step) = @_;
  if ($self->getStepByName($rO_step->getName())) {
    die "Error: step name \"" . $rO_step->getName() . "\" is not unique for this pipeline!";
  } else {
    push(@{$self->getSteps()}, $rO_step);
  }
}

sub getStepByName {
  my ($self, $name) = @_;

  my @namedSteps = grep($_->getName() eq $name, @{$self->getSteps()});

  # Should find 1 step at most since step names are unique
  if (@namedSteps) {
    return $namedSteps[0];
  } else {
    return undef;
  }
}

sub getSamples {
  my ($self) = @_;
  return $self->{'_samples'};
}

sub getNbSamples {
  my ($self) = @_;
  return scalar(@{$self->{'_samples'}});
}

sub addSample {
  my ($self, $rO_sample) = @_;
  if ($self->getSampleByName($rO_sample->getName())) {
    die "Error: sample name \"" . $rO_sample->getName() . "\" is not unique for this pipeline!";
  } else {
    push(@{$self->getSamples()}, $rO_sample);
  }
}

sub getSampleByName {
  my ($self, $name) = @_;

  my @namedSamples = grep($_->getName() eq $name, @{$self->getSamples()});

  # Should find 1 sample at most since sample names are unique
  if (@namedSamples) {
    return $namedSamples[0];
  } else {
    return undef;
  }
}

sub getRunTypes {
  my ($self) = @_;

  my $runTypes = [];

  foreach my $sample (@{$self->getSamples()}) {
    foreach my $readSet (@{$sample->getReadSets()}) {
      push(@$runTypes, $readSet->getRunType());
    }
  }

  # Remove duplicates and return unique runTypes
  return [keys(%{{map {$_ => 1}  @$runTypes}})];
}

1;
