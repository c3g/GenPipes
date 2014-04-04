#!/usr/bin/env perl

=head1 NAME

I<Step>

=head1 SYNOPSIS

Object used to hold information on a Step

=head1 DESCRIPTION

Object used to hold information on a Step


=head1 DEPENDENCY

=cut

package Step;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin";

# Dependencies
#-----------------------
use Data::Dumper;
use File::stat;
use List::Util qw(max);
use Job;
use LoadConfig;

# SUB
#-----------------------
sub new {
  my $class = shift;
  my $name = shift;
  my $loop = shift;
  my $rA_parentSteps = shift;

  # Step name is used in Bash $JOB_ID variable, hence only alphanumeric and "_" characters are allowed
  unless ($name =~ /^[a-zA-Z_]\w+$/) {
    die "Error: step name \"$name\" is invalid (should match [a-zA-Z_][a-zA-Z0-9_]+)!";
  }

  my $self = {
    '_name' => $name,
    '_loop' => $loop,
    '_parentSteps' => $rA_parentSteps,
    '_jobs' => [],
  };

  bless($self, $class);
  return $self;
}

sub getName {
  my ($self) = @_;
  return $self->{'_name'};
}

sub getLoop {
  my ($self) = @_;
  return $self->{'_loop'};
}

sub getParentSteps {
  my ($self) = @_;
  return $self->{'_parentSteps'};
}

sub getJobs {
  my ($self) = @_;
  return $self->{'_jobs'};
}

sub getNbJobs {
  my ($self) = @_;
  return scalar(@{$self->{'_jobs'}});
}

sub addJob {
  my ($self, $rO_job) = @_;
  push(@{$self->getJobs()}, $rO_job);
  $rO_job->setStep($self);
}

sub submitJob {
  my $self = shift;
  my $rH_cfg = shift;
  my $rO_job = shift;
  my $commandIdx = shift;

  unless ($rO_job->isUp2Date()) {
    $self->addJob($rO_job);

    # Job ID shell variables are suffixed by job count thus unique
    my $jobId = $self->getName() . "_JOB_" . $self->getNbJobs();
    $rO_job->setId($jobId);

    my $stepName = $self->getName();
    my $jobName = $rO_job->getName();


    unless ($commandIdx) {
      $commandIdx = 0
    };
    my $command = $rO_job->getCommand($commandIdx);

    # Concat dependencies job IDs with "$" as shell variables, and separated by cluster specific separator
    my $dependencies = join(LoadConfig::getParam($rH_cfg, $stepName, 'clusterDependencySep'), map("\$" . $_->getId(), @{$rO_job->getDependencies()}));

    # Print out job header and settings nicely
    my $separatorLine = "#" . "-" x 79 . "\n";
    print $separatorLine;
    print "# $jobId $jobName\n";
    print $separatorLine;
    print "JOB_NAME=$jobName\n";
    print "JOB_DEPENDENCIES=$dependencies\n";
    # Set job output filename based on job name and timestamp
    print "JOB_OUTPUT_RELATIVE_PATH=$stepName/\${JOB_NAME}_\$TIMESTAMP.o\n";
    print "JOB_OUTPUT=\$JOB_OUTPUT_ROOT/\$JOB_OUTPUT_RELATIVE_PATH\n";
    print "mkdir -p `dirname \$JOB_OUTPUT`\n";

    # Assign job number to job ID if any
    if (LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "true") {
      print $jobId . '=$(';
    }

    print "echo \"";

    my @doneFiles = map("$_.mugqic.done", @{$rO_job->getOutputFiles()});

    # Erase dones, on all jobs of the series
    if (@doneFiles > 0) {
      print "rm -f \\\n" . join(" \\\n", @doneFiles) . " && \\\n";
    }

    # Print out job modules
    if ($rO_job->getNbModules() > 0) {
      print "module load " . join(" ", @{$rO_job->getModules()}) . " && \\\n";
    }

    # Print out job command
    print "$command \\\n";
    print '&& echo \"MUGQICexitStatus:\$?\"';

    # Only add if it's the last job of the series.
    if (@doneFiles > 0 && $commandIdx == $rO_job->getNbCommands() - 1) {
      print " && touch \\\n" . join(" \\\n", @doneFiles);
    }
    print '"';
    print ' | ' . LoadConfig::getParam($rH_cfg, $stepName, 'clusterSubmitCmd');
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOtherArg');
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterWorkDirArg') . " \$WORK_DIR";
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOutputDirArg') . " \$JOB_OUTPUT";
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterJobNameArg') . " \$JOB_NAME";
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterWalltime');
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterQueue');
    print " " . LoadConfig::getParam($rH_cfg, $stepName, 'clusterCPU');
    if ($dependencies) {
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

    $rO_job->setCommandJobId($commandIdx, '$' . $jobId);

    # Write job parameters in job list file
    print "echo \"\$$jobId\t\$JOB_NAME\t\$JOB_DEPENDENCIES\t\$JOB_OUTPUT_RELATIVE_PATH\" >> \$JOB_LIST\n\n";
    return $jobId;
  }
}

1;
