#!/usr/env/perl

=head1 NAME

I<PipelineUtils>

=head1 SYNOPSIS

Utility classes for pipelines

=head1 DESCRIPTION

Utility classes for pipelines

=head1 DEPENDENCY

=cut

package PipelineUtils;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use LoadConfig;
use File::stat;
use Time::localtime;

# SUB
#-----------------------
sub testInputOutputs {
  my $rA_inputs = shift;
  my $rA_outputs = shift;

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
  my $runIt = 0;
  for my $outputFile (@{$rA_outputs}) {
    $retVal .= $outputFile.'.mugqic.done ';
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
    return $retVal;
  }
}

1;

