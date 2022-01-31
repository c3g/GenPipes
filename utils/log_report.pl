#!/usr/bin/env perl

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

=head1 NAME

I<Log>

=head1 SYNOPSIS

getLogReport <job_log_list_file>

=head1 DESCRIPTION

Parse internal log files and produce custom log reports

=head1 AUTHOR
B<Joel Fillon> - I<joel.fillon@mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

# Strict Pragmas
#---------------------
use strict;
use warnings;
#---------------------

# Dependencies
#-----------------------
use POSIX qw(strftime);
use File::Basename;
use File::Spec;
use Getopt::Long;

my $scriptName = basename($0);

my $usage=<<END;
NAME:
$scriptName

USAGE:
$scriptName [options] <job_list_file>

OPTIONS:
  -m,   --memtime       Output also memtime values if present in job output files
  -s,   --success       Show successful jobs only
  -nos, --nosuccess     Show unsuccessful jobs only i.e. failed or uncompleted jobs
  -h,   --help          Show this help

OUTPUT:
STDOUT

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Joel Fillon - joel.fillon\@mcgill.ca

END

# Check and assign the parameters
my ($memtimeOption, $successOption, $help);
GetOptions(
  "memtime" => \$memtimeOption,
  "success!" => \$successOption,
  "help"    => \$help
);
if ($help or not($ARGV[0])) {die $usage}

my $jobLogList = $ARGV[0];

print getLogTextReport($jobLogList);

# SUB
#-----------------------

# Return an array of hashes of log values for each job's log file created by the cluster
sub readJobLogListFile {
  my $jobLogListPath = shift;
  # Return array of hash of log values
  my $rAoH_jobLogList = shift;

  # Read the global list of log files
  open(JOB_LOG_LIST_FILE, $jobLogListPath) or die "Error: cannot open file $jobLogListPath\n" . $usage;
  while (my $line = <JOB_LOG_LIST_FILE>) {
    chomp($line);

    # Retrieve each job's log file path
    my ($jobId, $jobName, $jobDependencies, $clusterJobLogPath) = split(/\t/, $line);

    if (defined $jobName and defined $clusterJobLogPath) {
      my %jobLog;

      $jobLog{'jobId'} = $jobId;
      $jobLog{'jobName'} = $jobName;
      $jobLog{'jobDependencies'} = $jobDependencies;
      $jobLog{'memtime'} = ();

      # In old job log list version, cluster job log path was absolute. Now it is relative to job log list directory (useful if project directory has been moved or renamed).
      my $clusterJobLogFullPath = File::Spec->rel2abs(dirname($jobLogListPath)) . "/" . $clusterJobLogPath;
      unless (-e $clusterJobLogFullPath) {
        # Assume cluster job log path is absolute (old version)
        $clusterJobLogFullPath = $clusterJobLogPath;
      }
      $jobLog{'path'} = $clusterJobLogFullPath;

      # Read the job log file
      if (open(CLUSTER_JOB_LOG_FILE, $clusterJobLogFullPath)) {
        while (my $jobLine = <CLUSTER_JOB_LOG_FILE>) {
          # Job start date
          if ($jobLine =~ /^Begin PBS Prologue (.*) (\d+)$/) {
            $jobLog{'startDate'} = $1;
            $jobLog{'startSecondsSinceEpoch'} = $2;
          # Job number
          } elsif ($jobLine =~ /^Job ID:\s+(\S+)/) {
            $jobLog{'jobFullId'} = $1;
          # Job MUGQIC exit status
          } elsif ($jobLine =~ /MUGQICexitStatus:(\S+)/) {
            $jobLog{'MUGQICexitStatus'} = $1;
          # Username
          } elsif ($jobLine =~ /^Username:\s+(\S+)/) {
            $jobLog{'username'} = $1;
          # Group
          } elsif ($jobLine =~ /^Group:\s+(\S+)/) {
            $jobLog{'group'} = $1;
          # Session
          } elsif ($jobLine =~ /^Session:\s+(\S+)/) {
            $jobLog{'session'} = $1;
          # Limits
          } elsif ($jobLine =~ /^Limits:\s+(\S+)/) {
            $jobLog{'limits'} = $1;
          # Job used resources
          } elsif ($jobLine =~ /^Resources:\s+cput=(\S+),mem=(\S+),vmem=(\S+),walltime=(\d+:\d+:\d+)/) {
            $jobLog{'cput'} = $1;
            $jobLog{'mem'} = $2;
            $jobLog{'vmem'} = $3;
            $jobLog{'walltime'} = $4;
            # Compute duration in seconds from walltime hours, minutes, seconds
            $jobLog{'duration'} = timeToSeconds($4);
          # Queue
          } elsif ($jobLine =~ /^Queue:\s+(\S+)/) {
            $jobLog{'queue'} = $1;
          # Account
          } elsif ($jobLine =~ /^Account:\s+(\S+)/) {
            $jobLog{'account'} = $1;
          # Job exit status (should be the same as MUGQIC exit status unless MUGQIC exit status is skipped)
          # abacus syntax: "Exit_status", guillimin phase 2 syntax: "Exit code"
          } elsif ($jobLine =~ /^Exit(_status| code):\s+(\S+)/) {
            $jobLog{'exitStatus'} = $2;
          # Nodes
          } elsif ($jobLine =~ /^Nodes:\s+(\S+)/) {
            $jobLog{'nodes'} = $1;
          # Job end date
          } elsif ($jobLine =~ /^End PBS Epilogue (.*) (\d+)$/) {
            $jobLog{'endDate'} = $1;
            $jobLog{'endSecondsSinceEpoch'} = $2;
          } elsif ($jobLine =~ m/(\d+\.\d+) user, (\d+\.\d+) system, (\d+\.\d+) elapsed -- Max VSize = (\d+)KB, Max RSS = (\d+)KB/) {
            push @{$jobLog{'memtime'}}, {
              userTime => $1,
              systemTime => $2,
              elapsedTime => $3,
              maxVSize => $4,
              maxRSS => $5
            };
          }
        }
        close(CLUSTER_JOB_LOG_FILE);
      }

      # Filter jobs according to --success/--nosuccess option if any
      if (not(defined($successOption)) or    # Keep all jobs
          (defined($successOption) and
           # Keep successful jobs only
           ($successOption and exists($jobLog{'MUGQICexitStatus'}) and $jobLog{'MUGQICexitStatus'} == 0) or
           # Keep unsuccessful jobs only
           (not($successOption) and (not(exists($jobLog{'MUGQICexitStatus'})) or $jobLog{'MUGQICexitStatus'} != 0)))) {
        push (@$rAoH_jobLogList, \%jobLog);
      }
    }
  }
  close(JOB_LOG_LIST_FILE);
}

# Print out a log report in simple text format
sub getLogTextReport {

  my $jobLogListPath = shift;
  my @AoH_jobLogList;

  # Parse the job log files and get an array of hashes of those logs
  readJobLogListFile($jobLogListPath, \@AoH_jobLogList);

  my $logTextReport = "";

  $logTextReport .= "# Number of jobs: " . ($#AoH_jobLogList + 1) . "\n#\n";

  # Retrieve job status, first job start date, last job end date, shortest/longest jobs, lowest/highest memory jobs
  my $successfulJobCount = 0;
  my $activeJobCount = 0;
  my $inactiveJobCount = 0;
  my $failedJobCount = 0;
  my $firstStartSecondsSinceEpoch;
  my $lastEndSecondsSinceEpoch;
  my $shortestJob;
  my $longestJob;
  my $lowestMemoryJob;
  my $highestMemoryJob;

  for my $jobLog (@AoH_jobLogList) {
    if (exists $jobLog->{'MUGQICexitStatus'} and $jobLog->{'MUGQICexitStatus'} eq 0) {
      $successfulJobCount++;
      $jobLog->{'status'} = "SUCCESS";
    } elsif (exists $jobLog->{'endSecondsSinceEpoch'}) {
      # If job end date exists and MUGQICexitStatus != 0, job failed
      $failedJobCount++;
      $jobLog->{'status'} = "FAILED";
    } elsif (exists $jobLog->{'startSecondsSinceEpoch'}) {
      # If job start date exists but job end date does not, job is still running
      $activeJobCount++;
      $jobLog->{'status'} = "ACTIVE";
    } else {
      $inactiveJobCount++;
      $jobLog->{'status'} = "INACTIVE";
    }

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

    if (exists $jobLog->{'mem'}) {
      if (not defined $lowestMemoryJob or kiBToNum($lowestMemoryJob->{'mem'}) > kiBToNum($jobLog->{'mem'})) {
        $lowestMemoryJob = $jobLog;
      }
      if (not defined $highestMemoryJob or kiBToNum($highestMemoryJob->{'mem'}) < kiBToNum($jobLog->{'mem'})) {
        $highestMemoryJob = $jobLog;
      }
    }
  }

  # Print out job status counts
  $logTextReport .= "# Number of successful jobs: " . $successfulJobCount . "\n";
  $logTextReport .= "# Number of active jobs: " . $activeJobCount . "\n";
  $logTextReport .= "# Number of inactive jobs: " . $inactiveJobCount . "\n";
  $logTextReport .= "# Number of failed jobs: " . $failedJobCount . "\n#\n";

  # Print out execution time
  my $executionTime = (defined $firstStartSecondsSinceEpoch and defined $lastEndSecondsSinceEpoch) ? formatDuration($lastEndSecondsSinceEpoch - $firstStartSecondsSinceEpoch) : "N/A";
  my $startDate = defined $firstStartSecondsSinceEpoch ? strftime('%FT%T', localtime($firstStartSecondsSinceEpoch)) : "N/A";
  my $endDate = defined $lastEndSecondsSinceEpoch ? strftime('%FT%T', localtime($lastEndSecondsSinceEpoch)) : "N/A";
  $logTextReport .= "# Execution time: $startDate - $endDate ($executionTime)\n#\n";

  # Print out shortest and longest jobs
  $logTextReport .= "# Shortest job: " . (defined $shortestJob ? $shortestJob->{'jobName'} . " (" . formatDuration($shortestJob->{'duration'}) . ")" : "N/A") . "\n";
  $logTextReport .= "# Longest job: " . (defined $longestJob ? $longestJob->{'jobName'} . " (" . formatDuration($longestJob->{'duration'}) . ")" : "N/A") . "\n";
  $logTextReport .= "#\n";

  # Print out lowest and highest memory jobs
  $logTextReport .= "# Lowest memory job: " . (defined $lowestMemoryJob ? $lowestMemoryJob->{'jobName'} . " (" . sprintf("%.2f", kiBToGiB($lowestMemoryJob->{'mem'})) . " GiB" . ")" : "N/A") . "\n";
  $logTextReport .= "# Highest memory job: " . (defined $highestMemoryJob ? $highestMemoryJob->{'jobName'} . " (" . sprintf("%.2f", kiBToGiB($highestMemoryJob->{'mem'})) . " GiB" . ")" : "N/A") . "\n";
  $logTextReport .= "#\n";

  $logTextReport .= join("\t", (
    "#JOB_ID",
    "JOB_FULL_ID",
    "JOB_NAME",
    "JOB_DEPENDENCIES",
    "STATUS",
    "JOB_EXIT_CODE",
    "CMD_EXIT_CODE",
    "REAL_TIME",
    "START_DATE",
    "END_DATE",
    "CPU_TIME",
    "CPU_REAL_TIME_RATIO",
    "PHYSICAL_MEM",
    "VIRTUAL_MEM",
    "EXTRA_VIRTUAL_MEM_PCT",
    "LIMITS",
    "QUEUE",
    "USERNAME",
    "GROUP",
    "SESSION",
    "ACCOUNT",
    "NODES",
    "PATH"
  )) . "\n";

  if ($memtimeOption) {
    $logTextReport .= join("\t", (
      "#MEMTIME",
      "MEMTIME_USER",
      "MEMTIME_SYSTEM",
      "MEMTIME_ELAPSED",
      "MEMTIME_MAXVSIZE",
      "MEMTIME_MAXRSS"
    )) . "\n";
  }

  for my $jobLog (@AoH_jobLogList) {
    $logTextReport .= join("\t", (
      exists $jobLog->{'jobId'} ? $jobLog->{'jobId'} : "N/A",
      exists $jobLog->{'jobFullId'} ? $jobLog->{'jobFullId'} : "N/A",
      exists $jobLog->{'jobName'} ? $jobLog->{'jobName'} : "N/A",
      exists $jobLog->{'jobDependencies'} ? $jobLog->{'jobDependencies'} : "N/A",
      exists $jobLog->{'status'} ? $jobLog->{'status'} : "N/A",
      exists $jobLog->{'exitStatus'} ? $jobLog->{'exitStatus'} : "N/A",
      exists $jobLog->{'MUGQICexitStatus'} ? $jobLog->{'MUGQICexitStatus'} : "N/A",
      exists $jobLog->{'walltime'} ? $jobLog->{'walltime'} . " (" . formatDuration($jobLog->{'duration'}) . ")" : "N/A",
      exists $jobLog->{'startSecondsSinceEpoch'} ? strftime('%FT%T', localtime($jobLog->{'startSecondsSinceEpoch'})) : "N/A",
      exists $jobLog->{'endSecondsSinceEpoch'} ? strftime('%FT%T', localtime($jobLog->{'endSecondsSinceEpoch'})) : "N/A",
      exists $jobLog->{'cput'} ? $jobLog->{'cput'} . " (" . formatDuration(timeToSeconds($jobLog->{'cput'})) . ")" : "N/A",
      (exists $jobLog->{'walltime'} and exists $jobLog->{'cput'} and timeToSeconds($jobLog->{'walltime'}) != 0) ? sprintf("%.2f", timeToSeconds($jobLog->{'cput'}) / timeToSeconds($jobLog->{'walltime'})) : "N/A",
      exists $jobLog->{'mem'} ? sprintf("%.2f", kiBToGiB($jobLog->{'mem'})) . " GiB" : "N/A",
      exists $jobLog->{'vmem'} ? sprintf("%.2f", kiBToGiB($jobLog->{'vmem'})) . " GiB" : "N/A",
      (exists $jobLog->{'vmem'} and exists $jobLog->{'mem'} and kiBToGiB($jobLog->{'mem'}) != 0) ? sprintf("%.1f", (kiBToGiB($jobLog->{'vmem'}) / kiBToGiB($jobLog->{'mem'}) - 1) * 100) . " %" : "N/A",
      exists $jobLog->{'limits'} ? $jobLog->{'limits'} : "N/A",
      exists $jobLog->{'queue'} ? $jobLog->{'queue'} : "N/A",
      exists $jobLog->{'username'} ? $jobLog->{'username'} : "N/A",
      exists $jobLog->{'group'} ? $jobLog->{'group'} : "N/A",
      exists $jobLog->{'session'} ? $jobLog->{'session'} : "N/A",
      exists $jobLog->{'account'} ? $jobLog->{'account'} : "N/A",
      exists $jobLog->{'nodes'} ? $jobLog->{'nodes'} : "N/A",
      exists $jobLog->{'path'} ? $jobLog->{'path'} : "N/A"
    )) . "\n";

    # Memtime log on a separate line prefixed by "MEMTIME".
    if ($memtimeOption) {
      foreach my $memtime (@{$jobLog->{'memtime'}}) {
        $logTextReport .= join("\t", (
          "MEMTIME",
          exists $memtime->{'userTime'} ? formatDuration($memtime->{'userTime'}) : "N/A",
          exists $memtime->{'systemTime'} ? formatDuration($memtime->{'systemTime'}) : "N/A",
          exists $memtime->{'elapsedTime'} ? formatDuration($memtime->{'elapsedTime'}) : "N/A",
          exists $memtime->{'maxVSize'} ? sprintf("%.6f", kiBToGiB($memtime->{'maxVSize'})) . " GiB" : "N/A",
          exists $memtime->{'maxRSS'} ? sprintf("%.6f", kiBToGiB($memtime->{'maxRSS'})) . " GiB" : "N/A"
        )) . "\n";
      }
    }
  }
  return $logTextReport;
}

# Return seconds from time given in hh:mm:ss
sub timeToSeconds {
  my $time = shift;

  if ($time =~ /^(\d+):(\d\d):(\d\d)/) {
    return $1 * 60 ** 2 + $2 * 60 + $3;
  } else {
    return "N/A";
  }
}

# Return duration given in seconds into human readable format
sub formatDuration {
  my $seconds = shift;

  # Less than 1 minute
  if ($seconds < 60) {
    return $seconds . " s";
  }
  # Less than 1 hour
  elsif ($seconds < (60 * 60)) {
    return int($seconds / 60 ) . " min " . formatDuration($seconds % 60);
  }
  # Less than 1 day
  elsif ($seconds < (60 * 60 * 24)) {
    return int($seconds / (60 * 60)) . " h " . formatDuration($seconds % (60 * 60));
  }
  # 1 day or more
  else {
    return int($seconds / (60 * 60 * 24)) . " d " . formatDuration($seconds % (60 * 60 * 24));
  }
}

sub kiBToNum {
  my $size = shift;

  if ($size =~ /^(\d+)kb/i) {
    return $1;
  } else {
    return "N/A";
  }
}

sub kiBToGiB {
  my $size = shift;

  if ($size =~ /^(\d+)(kb)?/) {
    return $1 / (1024 ** 2);
  } else {
    return "N/A";
  }
}

1;
