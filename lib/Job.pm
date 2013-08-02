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
      }
      push(@{$self->{'_commands'}}, $command);
    }
    return $self->{'_commands'};
}

sub getCommands {
    my( $self ) = @_;
    return $self->{'_commands'};
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

1;
