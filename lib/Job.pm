#!/usr/env/perl

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

# Dependencies
#-----------------------
use File::stat;
use Time::localtime;
use LoadConfig;

# SUB
#-----------------------
sub new {
  my $class = shift;
  my $self = {
    '_isUp2Date' => 0,
  };
  bless($self, $class);
  return $self;
}

sub addCommand {
    my ( $self, $command ) = @_;
    if(defined($command)) {
      if(!defined($self->{'_commands'})) {
        $self->{'_commands'} = ();
        $self->{'_commandsJobId'} = ();
      }
      push(@{$self->{'_commands'}}, $command);
    }
}

sub addFilesToTest {
    my ( $self, $rA_filesToTest ) = @_;
    if(defined($rA_filesToTest)) {
      if(!defined($self->{'_filesToTest'})) {
        $self->{'_filesToTest'} = ();
      }
      push(@{$self->{'_filesToTest'}}, @{$rA_filesToTest});
    }
}

sub getFilesToTest {
  my ( $self ) = @_;
  return $self->{'_filesToTest'};
}

sub setOutputFileHash {
  my( $self, $rH_outputFiles ) = @_;
  $self->{'_outputFiles'} = $rH_outputFiles;
}

sub getOutputFileHash {
  my( $self ) = @_;
  return $self->{'_outputFiles'};
}

sub setCommandJobId {
  my( $self, $idx, $jobIdVarName ) = @_;
  $self->{'_commandsJobId'}->[$idx] = $jobIdVarName;
}

sub getCommandJobId {
  my( $self, $idx ) = @_;
  if(!defined($idx)) {
    $idx = 0;
  }
  return $self->{'_commandsJobId'}->[$idx];
}

sub getCommands {
    my( $self ) = @_;
    return $self->{'_commands'};
}

sub getNbCommands {
    my( $self ) = @_;
    return scalar(@{$self->{'_commands'}});
}

sub getCommand {
    my( $self, $idx ) = @_;
    if(!defined($idx)) {
      $idx = 0;
    }
    return $self->{'_commands'}->[$idx];
}

sub setUp2Date {
    my ( $self, $up2date ) = @_;
    if(defined($up2date)) {
      $self->{'_isUp2Date'} = $up2date;
    }
    return $self->{'_isUp2Date'};
}

sub isUp2Date {
    my( $self ) = @_;
    return $self->{'_isUp2Date'};
}

sub testInputOutputs {
  my( $self, $rA_inputs, $rA_outputs ) = @_;

  if(!defined($rA_inputs) || !defined($rA_outputs) || scalar(@{$rA_inputs}) == 0 || scalar(@{$rA_outputs}) == 0) {
    # Don't return touch, but return something so undef tests fail
    return "";
  }

  my $latestInput = ctime(stat($rA_inputs->[0])->mtime);
  for my $inputFile (@{$rA_inputs}) {
    my $inputTime = ctime(stat($inputFile)->mtime);
    if($inputTime > $latestInput) {
      $latestInput = $inputTime;
    }
  }

  my $retVal = " && touch ";
  my @filesToTest;
  my $runIt = 0;
  for my $outputFile (@{$rA_outputs}) {
    $retVal .= $outputFile.'.mugqic.done ';
    push(@filesToTest, $outputFile.'.mugqic.done');
    if($runIt != 0) {
      next;
    }

    if(!(-e $outputFile) || !(-e $outputFile.'.mugqic.done')) {
      $runIt = 1;
    }
    else {
      my $outputTime = ctime(stat($outputFile)->mtime);
      if($outputTime < $latestInput) {
        $runIt = 1;
      }
    }
  }

  if($runIt == 0) {
    return undef;
  }
  else {
    $self->setUp2Date(0);
    $self->addFilesToTest(\@filesToTest);
    return $retVal;
  }
}

1;
