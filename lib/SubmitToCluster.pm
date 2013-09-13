#!/usr/env/perl

=head1 NAME

I<SubmitToCluster>

=head1 SYNOPSIS

SubmitToCluster->printSubmitCmd()

=head1 DESCRIPTION

B<SubmitToCluster> is a library that reads from a config file and 
submits jobs to a cluster.

=head1 AUTHOR

B<Louis Letourneau> - I<louis.letourneau@mail.mcgill.ca>

=head1 DEPENDENCY

=cut

package SubmitToCluster;

# Strict Pragmas
#---------------------
use strict;
use warnings;

#---------------------

# Dependencies
#--------------------
use LoadConfig;

#--------------------

# SUB
#--------------------
sub initPipeline {
  #my $rH_cfg = shift;
  my $workDir = shift;

  # Set working directory to current one by default
  unless (defined $workDir and -d $workDir) {
    $workDir = "`pwd`";
  }
  print "cd $workDir\n";
  print "TIMESTAMP=`date +%FT%H.%M.%S`\n\n";
}

sub printSubmitCmd {
  my $rH_cfg = shift;
  my $stepName = shift;
  my $jobNameSuffix = shift;
  my $jobIdPrefix = shift;
  my $dependencyName = shift;
  my $sampleName = shift;
  my $command = shift;
  my $workDir = shift;

  # Set Job ID
  my $jobIdVarName = uc($jobIdPrefix) . "_JOB_ID";
  $jobIdVarName =~ s/\W/_/g;

  # Set working directory to current one by default
  unless (defined $workDir and -d $workDir) {
    $workDir = "`pwd`";
  }

  # Set job name and job log output directory depending on a global or sample-based step
  my $jobName = $stepName;
  my $jobOutputDir = $workDir . "/job_output/";
  if (defined($sampleName) and $sampleName ne "") {
    $jobName .= ".$sampleName";
    $jobOutputDir .= $sampleName;
  } else {
    $jobOutputDir .= "global";
  }
  if (defined($jobNameSuffix) and $jobNameSuffix ne "") {
    $jobName .= ".$jobNameSuffix";
  }

  # Set job log filename based on job name and timestamp
  my $jobOutputLog = "$jobOutputDir/$jobName" . "_\${TIMESTAMP}.o";

  print "mkdir -p $jobOutputDir\n";

  if (LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "true") {
    print $jobIdVarName . '=$(';
  }
  print 'echo "' . $command . ' && echo \"MUGQICexitStatus:\$?\"" | ';
  print LoadConfig::getParam($rH_cfg, $stepName, 'clusterSubmitCmd');
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOtherArg');
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterWorkDirArg') . " $workDir";
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOutputDirArg') . " $jobOutputLog";
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterJobNameArg') . " $jobName";
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterWalltime');
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterQueue');
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterCPU');
  if (defined($dependencyName)) {
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterDependencyArg') . $dependencyName;
  }
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterSubmitCmdSuffix');
  if (LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "true") {
    print ")";
  }
  print "\n";

  if (LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "false") {
    print "$jobIdVarName=$jobName\n";
  }

  # Write job ID, name and log path in job list file
  print "echo \"\${$jobIdVarName}\t$jobName\t$jobOutputLog\" >> $workDir/job_output/job_list_\${TIMESTAMP}\n\n";

  return $jobIdVarName;
}
1;
