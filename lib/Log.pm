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
use POSIX qw(strftime);

use LoadConfig;

# SUB
#-----------------------

# Return a hash indexed by job ID of hash of log values for each job's log file created by the cluster
sub readJobLogListFile {
  my $rH_cfg = shift;
  my $workDirectory = shift;
  my $jobLogListPath = shift;
  # Return hash of hash of log values
  my $rHoH_jobLogList = shift;

  # Read the global list of log files
  open(JOB_LOG_LIST_FILE, $jobLogListPath) or die "Cannot open $jobLogListPath\n";
  while(my $line = <JOB_LOG_LIST_FILE>) {
    # Retrieve each job's log file path
    if($line =~ /^(\S+);(\S+)/) {
      my $jobId = $1;
      my $clusterJobLogPath = $workDirectory . $2;
  
      # Read the job log file
      if (open(CLUSTER_JOB_LOG_FILE, $clusterJobLogPath)) {
        $rHoH_jobLogList->{$jobId}{'path'} = $clusterJobLogPath;
        while(my $jobLine = <CLUSTER_JOB_LOG_FILE>) {
          # Job start date
          if($jobLine =~ /^Begin PBS Prologue (.*) (\d+)$/) {
            $rHoH_jobLogList->{$jobId}{'startDate'} = $1;
            $rHoH_jobLogList->{$jobId}{'startSecondsSinceEpoch'} = $2;
          # Job number
          } elsif($jobLine =~ /^Job ID:\s+(\S+)/) {
            $rHoH_jobLogList->{$jobId}{'jobNumber'} = $1;
          # Job name
          } elsif($jobLine =~ /^Job Name:\s+(\S+)/) {
            $rHoH_jobLogList->{$jobId}{'jobName'} = $1;
          # Job exit status
          #} elsif($jobLine =~ /^MUGQICexitStatus:(\d+)/) {
          } elsif($jobLine =~ /^Exit_status:\s+(\d+)/) {
            $rHoH_jobLogList->{$jobId}{'MUGQICexitStatus'} = $1;
          # Job used resources
          } elsif($jobLine =~ /^Resources:\s+cput=(\S+),mem=(\S+),vmem=(\S+),walltime=((\d+):(\d+):(\d+))/) {
            $rHoH_jobLogList->{$jobId}{'cput'} = $1;
            $rHoH_jobLogList->{$jobId}{'mem'} = $2;
            $rHoH_jobLogList->{$jobId}{'vmem'} = $3;
            $rHoH_jobLogList->{$jobId}{'walltime'} = $4;
            $rHoH_jobLogList->{$jobId}{'duration'} = $5 * 60 * 60 + $6 * 60 + $7;
          # Job end date
          } elsif($jobLine =~ /^End PBS Epilogue (.*) (\d+)$/) {
            $rHoH_jobLogList->{$jobId}{'endDate'} = $1;
            $rHoH_jobLogList->{$jobId}{'endSecondsSinceEpoch'} = $2;
          }
        }
        close(CLUSTER_JOB_LOG_FILE);
      } else {
        warn "Cannot open $clusterJobLogPath\n";
      }
    }
  }
  close(JOB_LOG_LIST_FILE);
}

sub getLogTextReport {
  my $rHoH_jobLogList = shift;
  my $logTextReport = "";

  # Retrieve start and dates for the whole job pipeline
  my @A_jobSortByStartDate = sort {$rHoH_jobLogList->{$a}{'startSecondsSinceEpoch'} <=> $rHoH_jobLogList->{$b}{'startSecondsSinceEpoch'}} keys %$rHoH_jobLogList;
  my @A_jobSortByEndDate = sort {$rHoH_jobLogList->{$a}{'endSecondsSinceEpoch'} <=> $rHoH_jobLogList->{$b}{'endSecondsSinceEpoch'}} keys %$rHoH_jobLogList;

  my %startJob = %{$rHoH_jobLogList->{$A_jobSortByStartDate[0]}};
  my $startSecondsSinceEpoch = $startJob{'startSecondsSinceEpoch'};
  my %endJob = %{$rHoH_jobLogList->{$A_jobSortByEndDate[$#A_jobSortByEndDate]}};
  my $endSecondsSinceEpoch = $endJob{'endSecondsSinceEpoch'};
  $logTextReport .= "Execution time: " . strftime('%FT%T', localtime($startSecondsSinceEpoch)) . " - " . strftime('%FT%T', localtime($endSecondsSinceEpoch)) . " (" . formatDuration($endSecondsSinceEpoch - $startSecondsSinceEpoch) . ")\n\n";

  # Retrieve shortest and longest jobs
  my @A_jobSortByDuration = sort {$rHoH_jobLogList->{$a}{'duration'} <=> $rHoH_jobLogList->{$b}{'duration'}} keys %$rHoH_jobLogList;
  my %shortestJob = %{$rHoH_jobLogList->{$A_jobSortByDuration[0]}};
  my %longestJob = %{$rHoH_jobLogList->{$A_jobSortByDuration[$#A_jobSortByDuration]}};
  $logTextReport .= "Shortest job: " . $shortestJob{'jobName'} . " (" . $shortestJob{'walltime'} . ")\n";
  $logTextReport .= "Longest job: " . $longestJob{'jobName'} . " (" . $longestJob{'walltime'} . ")\n\n";

  $logTextReport .= join("\t", (
    "JOB_NAME",
    "EXIT_CODE",
    "WALL_TIME",
    "START_DATE",
    "END_DATE",
    "CPU_TIME",
    "MEMORY",
    "VMEMORY",
  )) . "\n";

  for my $jobLog (sort {$rHoH_jobLogList->{$a}{'startSecondsSinceEpoch'} <=> $rHoH_jobLogList->{$b}{'startSecondsSinceEpoch'}} keys %$rHoH_jobLogList) {
    $logTextReport .= join("\t", (
      $rHoH_jobLogList->{$jobLog}{'jobName'},
      $rHoH_jobLogList->{$jobLog}{'MUGQICexitStatus'},
      $rHoH_jobLogList->{$jobLog}{'walltime'},
      strftime('%FT%T', localtime($rHoH_jobLogList->{$jobLog}{'startSecondsSinceEpoch'})),
      strftime('%FT%T', localtime($rHoH_jobLogList->{$jobLog}{'endSecondsSinceEpoch'})),
      $rHoH_jobLogList->{$jobLog}{'cput'},
      $rHoH_jobLogList->{$jobLog}{'mem'},
      $rHoH_jobLogList->{$jobLog}{'vmem'}
    )) . "\n";
  }

  return $logTextReport;
}

sub formatDuration
{
  my $seconds =  shift;

  if ($seconds < 60) {
    # less than a minute
    return $seconds . " s";
  }
  elsif ($seconds <= (60 * 60)) {
    # less than an hour
    return int($seconds / 60 ) . " min " . formatDuration($seconds % 60);
  }
  elsif ($seconds <= (60 * 60 * 24)) {
    # less than a day
    return int($seconds / (60 * 60)) . " h " . formatDuration($seconds % (60 * 60));
  }
  elsif ($seconds <= (60 * 60 * 24 * 7)) {
    # less than a week
    return int($seconds / (60 * 60 * 24)) . " day(s) " . formatDuration($seconds % (60 * 60 * 24));
  }
  else {
    # fall-back weeks ago
    return int($seconds / (60 * 60 * 24 * 7)) . " week(s) " . formatDuration($seconds % (60 * 60 * 24 * 7));
  }
}
        
1;
