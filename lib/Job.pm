#!/usr/bin/env perl

=head1 NAME

I<Job>

=head1 SYNOPSIS

Object used to hold information on a Job or Jobs to run

=head1 DESCRIPTION

Object used to hold information on a Job or Jobs to run


=head1 DEPENDENCY

=cut

package Job;

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
use Data::Dumper;
use File::stat;
use List::Util qw(max);
use LoadConfig;

# SUB
#-----------------------
sub new {
  my $class = shift;
  my $rA_inputFiles = shift;
  my $rA_outputFiles = shift;

  my $self = {
    '_inputFiles' => [],
    '_outputFiles' => []
  };
  bless($self, $class);

  foreach my $inputFile (@$rA_inputFiles) {
    if ($inputFile) {
      push (@{$self->getInputFiles()}, $inputFile);
    }
  }

  foreach my $outputFile (@$rA_outputFiles) {
    if ($outputFile) {
      push (@{$self->getOutputFiles()}, $outputFile);
    }
  }

  return $self;
}

sub getInputFiles {
  my ($self) = @_;
  return $self->{'_inputFiles'};
}

sub getNbInputFiles {
  my ($self) = @_;
  return scalar(@{$self->{'_inputFiles'}});
}

sub getOutputFiles {
  my ($self) = @_;
  return $self->{'_outputFiles'};
}

sub getNbOutputFiles {
  my ($self) = @_;
  return scalar(@{$self->{'_outputFiles'}});
}

# Job name reflect job step/loopTags hierarchy e.g <step>.<sampleLoopTag>.<readSetLoopTag>
sub getName {
  my ($self) = @_;
  return join(".", ($self->getStep()->getName(), @{$self->getLoopTags()}));
}

sub setId {
  my ($self, $id) = @_;
  $self->{'_id'} = $id;
}

sub getId {
  my ($self) = @_;
  return $self->{'_id'};
}

sub addCommand {
  my ($self, $command) = @_;
  if (defined($command)) {
    if (!defined($self->{'_commands'})) {
      $self->{'_commands'} = ();
      $self->{'_commandsJobId'} = ();
    }
    push(@{$self->{'_commands'}}, $command);
  }
}

sub addFilesToTest {
  my ($self, $rA_filesToTest) = @_;
  if (defined($rA_filesToTest)) {
    if (!defined($self->{'_filesToTest'})) {
      $self->{'_filesToTest'} = ();
    }
    push(@{$self->{'_filesToTest'}}, @{$rA_filesToTest});
  }
}

sub getFilesToTest {
  my ($self) = @_;
  return $self->{'_filesToTest'};
}

sub setOutputFileHash {
  my ($self, $rH_outputFiles) = @_;
  $self->{'_outputFiles'} = $rH_outputFiles;
}

sub getOutputFileHash {
  my ($self) = @_;
  return $self->{'_outputFiles'};
}

sub setCommandJobId {
  my ($self, $idx, $jobIdVarName) = @_;
  $self->{'_commandsJobId'}->[$idx] = $jobIdVarName;
}

sub getCommandJobId {
  my ($self, $idx) = @_;
  if (!defined($idx)) {
    $idx = 0;
  }
  return $self->{'_commandsJobId'}->[$idx];
}

sub getCommands {
  my ($self) = @_;
  return $self->{'_commands'};
}

sub getNbCommands {
  my ($self) = @_;
  return scalar(@{$self->{'_commands'}});
}

sub getCommand {
  my ($self, $idx) = @_;
  if (!defined($idx)) {
    $idx = 0;
  }
  return $self->{'_commands'}->[$idx];
}

sub setStep {
  my ($self, $rH_step) = @_;
  $self->{'_step'} = $rH_step;
}

sub getStep {
  my ($self) = @_;
  return $self->{'_step'};
}

sub setLoopTags {
  my ($self, $rA_loopTags) = @_;
  $self->{'_loopTags'} = $rA_loopTags;
}

sub getLoopTags {
  my ($self) = @_;
  return $self->{'_loopTags'};
}

# Test if job is up to date by checking job output .done presence and comparing input/output file modification times 
# If job is out of date, set job files to test with list of output files
sub isUp2Date {
  my ($self) = @_;

  if ($self->getNbInputFiles() == 0 || $self->getNbOutputFiles() == 0) {
    # Don't return 'touch' command, but return something so undef tests fail
    return "";
  }

  my $isJobUpToDate = 1;   # Job is up to date by default

  # Retrieve latest input file modification time i.e. maximum stat mtime
  # Use 'echo' system command to expand environment variables in input file paths if any
  # Also check if input file exists before calling mtime function, return 0 otherwise
  my $latestInputTime = max(map(-e `echo -n $_` ? stat(`echo -n $_`)->mtime : 0, @{$self->getInputFiles()}));
#warn "Job.isUp2Date: after latestInputTime $latestInputTime isJobUpToDate $isJobUpToDate\n";

  if ($latestInputTime == 0) {    # i.e. if job input files don't exist yet
    $isJobUpToDate = 0;
  } else {
    foreach my $outputFile (@{$self->getOutputFiles()}) {

#warn "Job.isUp2Date: outputFile $outputFile\n";
      # Use 'echo' system command to expand environment variables in output file path if any
      my $outputExpandedFile = `echo -n $outputFile`;
  
      # Skip further tests if job is already out of date
      if ($isJobUpToDate) {
        # If .done file is missing or if output file is older than latest input file, job is not up to date
        unless ((-e $outputExpandedFile) and (-e $outputExpandedFile . ".mugqic.done") and (stat($outputExpandedFile)->mtime >= $latestInputTime)) {
#warn "Job.isUp2Date: in unless outputExpandedFile $outputExpandedFile mtime " . stat($outputExpandedFile)->mtime . "\n";
          $isJobUpToDate = 0;
        }
      }
    }
  }

#warn "Job.isUp2Date: isJobUpToDate $isJobUpToDate\n";
  return $isJobUpToDate;

#  if ($isJobUpToDate) {
#    return undef;
#  } else {
#    my @outputDoneFiles = map("$_.mugqic.done", @$rA_outputs);
#    $self->addFilesToTest(\@outputDoneFiles);
#    return " && touch " . join(" ", @outputDoneFiles);
#  }
}

sub getDependencies {
  my ($self) = @_;

  my $dependencyJobs = [];

  foreach my $parentStep (@{$self->getStep()->getParentSteps()}) {
    foreach my $parentJob (@{$parentStep->getJobs()}) {
      my $isParentJobDependency = 1;    # By default parent job is a dependency
      
      # Compare job loop tags: parent job is a dependency if current job loop tags are
      # identical to parent job loop tags one by one or if loop tags are missing
      # on any side e.g:

      # current job loop tags [sampleA, readset2] are dependent of
      # parent job loop tags [sampleA, readset2] or [sampleA] or []
      # but not dependent of [sampleB, readset1] nor [sampleB].

      # Current job loop tags [sampleA] are dependent of
      # parent job loop tags [sampleA, readset2] or [sampleA] or []
      # but not [sampleB].

      # Current job empty loop tags [] are dependent of
      # every parent job loop tags [sampleA, readset2], [sampleA], [sample2], []...

      my @currentJobLoopTags = @{$self->getLoopTags()};
      my @parentJobLoopTags = @{$parentJob->getLoopTags()};
      foreach my $i (0 .. $#currentJobLoopTags) {
        if ($parentJobLoopTags[$i] and $currentJobLoopTags[$i] ne $parentJobLoopTags[$i]) {
          $isParentJobDependency = 0;
        }
      }

      if ($isParentJobDependency) {
        push(@$dependencyJobs, $parentJob);
      }
    }
  }

  return $dependencyJobs;
}

1;
