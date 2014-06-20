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
use lib $FindBin::Bin;

# Dependencies
#-----------------------
use Digest::MD5 qw(md5_hex);
use File::Basename;
use File::stat;
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
    if($rA_filesToTest) {
      if(!defined($self->{'_filesToTest'})) {
        $self->{'_filesToTest'} = ();
      }
      # Create checksum file name of all files to test and store it in the directory of the first file to test
      my $checksumFileToTest = dirname(@{$rA_filesToTest}[0]) . "/" . md5_hex(join("", @{$rA_filesToTest})) . ".mugqic.done";
      push(@{$self->{'_filesToTest'}}, $checksumFileToTest);
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

  $self->addFilesToTest($rA_outputs);

  my $rA_inputFiles = [];
  my $rA_outputFiles = [];

  # Update file paths with expanded environment variables if any
  for my $filePath (@$rA_inputs) {
    my $longFilePath = `echo $filePath`;
    chomp($longFilePath);
    push($rA_inputFiles, $longFilePath);
  }

  for my $filePath (@$rA_outputs) {
    my $longFilePath = `echo $filePath`;
    chomp($longFilePath);
    push($rA_outputFiles, $longFilePath);
  }

  my $latestInput = -1;
  my $inputStat = stat($rA_inputFiles->[0]);
  if(defined($inputStat)) {
    $latestInput = $inputStat->mtime;
  }
  for my $inputFile (@{$rA_inputFiles}) {
    $inputStat = stat($inputFile);
    my $inputTime = -1;
    if(defined($inputStat)) {
      $inputTime = $inputStat->mtime;
    }
    if($inputTime > $latestInput) {
      $latestInput = $inputTime;
    }
  }

  my $retVal = " && touch ";
  my $runIt = 0;

  for my $fileToTest (@{$self->getFilesToTest()}) {
    if(!(-e $fileToTest)) {
      $runIt = 1;
    }
  }

  for my $outputFile (@{$rA_outputFiles}) {
    $retVal .= $outputFile.'.mugqic.done ';
    if($runIt != 0) {
      next;
    }

    if(!(-e $outputFile)) {
      $runIt = 1;
    }
    else {
      my $outputTime = $latestInput-1; # make it older in case it doesn't exist
      my $outputStat = stat($outputFile);
      if(defined($outputStat)) {
        $outputTime = $outputStat->mtime;
      }
      if($outputTime < $latestInput) {
        $runIt = 1;
      }
    }
  }

  $self->setUp2Date(!$runIt);
  if($runIt == 0) {
    return undef;
  }
  else {
    return $retVal;
  }
}

1;
