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
  my $workDir = shift;

  # Set working directory to current one by default
  unless (defined $workDir and -d $workDir) {
    $workDir = "`pwd`";
  }

  # Set pipeline header and global variables
  print "#!/bin/bash\n\n";
  print "WORK_DIR=$workDir\n";
  print "JOB_OUTPUT_ROOT=\$WORK_DIR/job_output\n";
  print "TIMESTAMP=`date +%FT%H.%M.%S`\n";
  print "JOB_LIST=\$JOB_OUTPUT_ROOT/job_list_\$TIMESTAMP\n\n";
  print "cd \$WORK_DIR\n\n";
}

sub printSubmitCmd {
  my $rH_cfg = shift;
  my $stepName = shift;
  my $jobNameSuffix = shift;
  my $jobIdPrefix = shift;
  my $dependencies = shift;
  my $sampleName = shift;
  my $rO_job = shift;
  my $commandIdx = shift;

  if($rO_job->isUp2Date()) {
    return undef;
  }

  if(!defined($commandIdx)) {
    $commandIdx = 0;
  }
  $command = $rO_job->getCommand($commandIdx);

  # Set Job ID
  my $jobId = uc($jobIdPrefix) . "_JOB_ID";
  $jobId =~ s/\W/_/g;

  # Set job name and job output directory depending on a global or sample-based step
  my $jobName = $stepName;
  my $jobOutputDir = "\$JOB_OUTPUT_ROOT/";
  if (defined($sampleName) and $sampleName ne "") {
    $jobName .= ".$sampleName";
    $jobOutputDir .= $sampleName;
  } else {
    $jobOutputDir .= "global";
  }
  if (defined($jobNameSuffix) and $jobNameSuffix ne "") {
    $jobName .= ".$jobNameSuffix";
  }

  # Check if $dependencies is initialized
  unless (defined($dependencies)) {
    $dependencies = "";
  }

  # Print out job header and settings nicely
  print "#--------------------------------------------------------------------------------\n";
  print "# $jobId $jobName\n";
  print "#--------------------------------------------------------------------------------\n";
  print "JOB_NAME=$jobName\n";
  print "JOB_DEPENDENCIES=$dependencies\n";
  print "JOB_OUTPUT_DIR=$jobOutputDir\n";
  # Set job output filename based on job name and timestamp
  print "JOB_OUTPUT=\$JOB_OUTPUT_DIR/\${JOB_NAME}_\$TIMESTAMP.o\n";
  print "mkdir -p \$JOB_OUTPUT_DIR\n";

  # Assign job number to job ID if any
  if (LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "true") {
    print $jobId . '=$(';
  }
  # Print out job command
  print 'echo "' . $command . ' && echo \"MUGQICexitStatus:\$?\"" | ';
  print LoadConfig::getParam($rH_cfg, $stepName, 'clusterSubmitCmd');
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOtherArg');
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterWorkDirArg') . " \$WORK_DIR";
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOutputDirArg') . " \$JOB_OUTPUT";
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterJobNameArg') . " \$JOB_NAME";
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterWalltime');
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterQueue');
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterCPU');
  if ($dependencies ne "") {
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterDependencyArg') . "\$JOB_DEPENDENCIES";
  }
  print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterSubmitCmdSuffix');
  if (LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "true") {
    print ")";
  }
  print "\n";

  if (LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "false") {
    print "$jobId=$jobName\n";
  }

  $rO_job->setCommandJobId($commandIdx, '$'.$jobId);

  # Write job parameters in job list file
  print "echo \"\$$jobId\t\$JOB_NAME\t\$JOB_DEPENDENCIES\t\$JOB_OUTPUT\" >> \$JOB_LIST\n\n"; 
  return $jobId;
}
1;
