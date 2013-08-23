#!/usr/env/perl

=head1 NAME

I<Log>

=head1 SYNOPSIS

Log-> readRestartFile()

=head1 DESCRIPTION

B<Log> is a library to parse internal log files and produce custom log reports



=head1 AUTHOR
B<Joel Fillon> - I<TO_BE_ADDED>
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package Log;

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

# Return a hash indexed by job ID of hash of log values for each job's log file created by the cluster
sub readRestartFile {
  my $rH_cfg = shift;
  my $workDirectory = shift;
  my $restartLogFilePath = shift;
  # Return hash of hash of log values
  my $rHoH_restartList = shift;

  # Read the global restart log file
  open(RESTART_LOG_FILE, $restartLogFilePath) or die "Cannot open ".$restartLogFilePath."\n";
  while(my $line = <RESTART_LOG_FILE>) {
    # Retrieve each job's log file path
    if($line =~ /^(\S+);(\S+)/) {
      my $jobId = $1;
      my $mugqicJobLogFilePath = $workDirectory.$2;

      $rHoH_restartList->{$jobId}->{'path'} = $mugqicJobLogFilePath;

      # Read the job log file
      open(MUGQIC_JOB_LOG_FILE, $mugqicJobLogFilePath) or die "Cannot open ".$mugqicJobLogFilePath."\n"; 
      while(my $jobLine = <MUGQIC_JOB_LOG_FILE>) {
        # Job start date
        if($jobLine =~ /^Begin PBS Prologue (.*) \d+$/) {
          $rHoH_restartList->{$jobId}->{'startDate'} = $1;
        # Job number
        } elsif($jobLine =~ /^Job ID:\s+(\S+)/) {
          $rHoH_restartList->{$jobId}->{'jobNumber'} = $1;
        # Job name
        } elsif($jobLine =~ /^Job Name:\s+(\S+)/) {
          $rHoH_restartList->{$jobId}->{'jobName'} = $1;
        # Job exit status
        } elsif($jobLine =~ /^MUGQICexitStatus:(\d+)/) {
          $rHoH_restartList->{$jobId}->{'MUGQICexitStatus'} = $1;
        # Job used resources
        } elsif($jobLine =~ /^Resources:\s+cput=(\S+),mem=(\S+),vmem=(\S+),walltime=(\S+)/) {
          $rHoH_restartList->{$jobId}->{'cput'} = $1;
          $rHoH_restartList->{$jobId}->{'mem'} = $2;
          $rHoH_restartList->{$jobId}->{'vmem'} = $3;
          $rHoH_restartList->{$jobId}->{'walltime'} = $4;
        # Job end date
        } elsif($jobLine =~ /^End PBS Epilogue (.*) \d+$/) {
          $rHoH_restartList->{$jobId}->{'endDate'} = $1;
        }
      }
      close(MUGQIC_JOB_LOG_FILE);
    }
  }
  close(RESTART_LOG_FILE);

}

1;
