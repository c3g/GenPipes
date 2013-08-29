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

# Return an array of hashes of log values for each job's log file created by the cluster
sub readJobLogListFile {
  my $rH_cfg = shift;
  my $workDirectory = shift;
  my $jobLogListPath = shift;
  # Return array of hash of log values
  my $rAoH_jobLogList = shift;

  # Read the global list of log files
  open(JOB_LOG_LIST_FILE, $jobLogListPath) or die "Cannot open $jobLogListPath\n";
  while(my $line = <JOB_LOG_LIST_FILE>) {
    # Retrieve each job's log file path
    if($line =~ /^(\S+);(\S+)/) {
      my $jobId = $1;
      my $clusterJobLogPath = $workDirectory . $2;
      my %jobLog;
  
      $jobLog{'jobId'} = $jobId;
      $jobLog{'path'} = $clusterJobLogPath;
      # Read the job log file
      if (open(CLUSTER_JOB_LOG_FILE, $clusterJobLogPath)) {
        while(my $jobLine = <CLUSTER_JOB_LOG_FILE>) {
          # Job start date
          if($jobLine =~ /^Begin PBS Prologue (.*) (\d+)$/) {
            $jobLog{'startDate'} = $1;
            $jobLog{'startSecondsSinceEpoch'} = $2;
          # Job number
          } elsif($jobLine =~ /^Job ID:\s+(\S+)/) {
            $jobLog{'jobNumber'} = $1;
          # Job name
          } elsif($jobLine =~ /^Job Name:\s+(\S+)/) {
            $jobLog{'jobName'} = $1;
          # Job MUGQIC exit status
          } elsif($jobLine =~ /^MUGQICexitStatus:(\d+)/) {
            $jobLog{'MUGQICexitStatus'} = $1;
          # Job exit status (should be the same as MUGQIC exit status unless MUGQIC exit status is skipped)
          } elsif($jobLine =~ /^Exit_status:\s+(\d+)/) {
            $jobLog{'exitStatus'} = $1;
          # Job used resources
          } elsif($jobLine =~ /^Resources:\s+cput=(\S+),mem=(\S+),vmem=(\S+),walltime=((\d+):(\d+):(\d+))/) {
            $jobLog{'cput'} = $1;
            $jobLog{'mem'} = $2;
            $jobLog{'vmem'} = $3;
            $jobLog{'walltime'} = $4;
            # Compute duration in seconds from walltime hours, minutes, seconds
            $jobLog{'duration'} = $5 * 60 * 60 + $6 * 60 + $7;
          # Job end date
          } elsif($jobLine =~ /^End PBS Epilogue (.*) (\d+)$/) {
            $jobLog{'endDate'} = $1;
            $jobLog{'endSecondsSinceEpoch'} = $2;
          }
        }
        close(CLUSTER_JOB_LOG_FILE);
      }
      push (@$rAoH_jobLogList, \%jobLog);
    }
  }
  close(JOB_LOG_LIST_FILE);
}

sub getLogTextReport {
  my $rAoH_jobLogList = shift;
  my $logTextReport = "";

  # Retrieve first job start date, last job end date, shortest job, longest job
  my $firstStartSecondsSinceEpoch;
  my $lastEndSecondsSinceEpoch;
  my $shortestJob;
  my $longestJob;

  for my $jobLog (@$rAoH_jobLogList) {
    if (exists $jobLog->{'startSecondsSinceEpoch'} and (not defined $firstStartSecondsSinceEpoch or $firstStartSecondsSinceEpoch > $jobLog->{'startSecondsSinceEpoch'})) {
      $firstStartSecondsSinceEpoch = $jobLog->{'startSecondsSinceEpoch'};
    }
    if (exists $jobLog->{'endSecondsSinceEpoch'} and (not defined $lastEndSecondsSinceEpoch or $lastEndSecondsSinceEpoch < $jobLog->{'endSecondsSinceEpoch'})) {
      $lastEndSecondsSinceEpoch = $jobLog->{'endSecondsSinceEpoch'};
    }
    if (exists $jobLog->{'duration'}) {
      if (not defined $shortestJob or $shortestJob->{'duration'} > $jobLog->{'duration'}) {
        $shortestJob = $jobLog;
      }
      if (not defined $longestJob or $longestJob->{'duration'} < $jobLog->{'duration'}) {
        $longestJob = $jobLog;
      }
    }
  }

  # Print out execution time
  my $executionTime = (defined $firstStartSecondsSinceEpoch and defined $lastEndSecondsSinceEpoch) ? formatDuration($lastEndSecondsSinceEpoch - $firstStartSecondsSinceEpoch) : "N/A";
  my $startDate = defined $firstStartSecondsSinceEpoch ? strftime('%FT%T', localtime($firstStartSecondsSinceEpoch)) : "N/A";
  my $endDate = defined $lastEndSecondsSinceEpoch ? strftime('%FT%T', localtime($lastEndSecondsSinceEpoch)) : "N/A";
  $logTextReport .= "Execution time: $startDate - $endDate ($executionTime)\n\n";

  # Print out shortest and longest jobs
  $logTextReport .= "Shortest job: " . (defined $shortestJob ? $shortestJob->{'jobName'} . " (" . formatDuration($shortestJob->{'duration'}) . ")" : "N/A") . "\n";
  $logTextReport .= "Longest job: " . (defined $longestJob ? $longestJob->{'jobName'} . " (" . formatDuration($longestJob->{'duration'}) . ")" : "N/A") . "\n\n";

  $logTextReport .= join("\t", (
    "JOB_ID",
    "JOB_NAME",
    "EXIT_CODE",
    "WALL_TIME",
    "START_DATE",
    "END_DATE",
    "CPU_TIME",
    "MEMORY",
    "VMEMORY"
  )) . "\n";

  for my $jobLog (@$rAoH_jobLogList) {
    $logTextReport .= join("\t", (
      exists $jobLog->{'jobId'} ? $jobLog->{'jobId'} : "N/A",
      exists $jobLog->{'jobName'} ? $jobLog->{'jobName'} : "N/A",
      exists $jobLog->{'exitStatus'} ? $jobLog->{'exitStatus'} : "N/A",
      exists $jobLog->{'walltime'} ? $jobLog->{'walltime'} : "N/A",
      exists $jobLog->{'startSecondsSinceEpoch'} ? strftime('%FT%T', localtime($jobLog->{'startSecondsSinceEpoch'})) : "N/A",
      exists $jobLog->{'endSecondsSinceEpoch'} ? strftime('%FT%T', localtime($jobLog->{'endSecondsSinceEpoch'})) : "N/A",
      exists $jobLog->{'cput'} ? $jobLog->{'cput'} : "N/A",
      exists $jobLog->{'mem'} ? $jobLog->{'mem'} : "N/A",
      exists $jobLog->{'vmem'} ? $jobLog->{'vmem'} : "N/A"
    )) . "\n";
  }

  return $logTextReport;
}

# Return duration given in seconds into human readable format
sub formatDuration {
  my $seconds =  shift;

  # Less than 1 minute
  if ($seconds < 60) {
    return $seconds . " s";
  }
  # Less than 1 hour
  elsif ($seconds <= (60 * 60)) {
    return int($seconds / 60 ) . " min " . formatDuration($seconds % 60);
  }
  # Less than 1 day
  elsif ($seconds <= (60 * 60 * 24)) {
    return int($seconds / (60 * 60)) . " h " . formatDuration($seconds % (60 * 60));
  }
  # Less than 1 week
  elsif ($seconds <= (60 * 60 * 24 * 7)) {
    return int($seconds / (60 * 60 * 24)) . " day(s) " . formatDuration($seconds % (60 * 60 * 24));
  }
  # 1 week or more
  else {
    return int($seconds / (60 * 60 * 24 * 7)) . " week(s) " . formatDuration($seconds % (60 * 60 * 24 * 7));
  }
}
        
1;
