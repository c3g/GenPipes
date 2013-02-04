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
sub initSubmit {
    my $rH_cfg     = shift;
    my $sampleName = shift;

    print "mkdir -p " . LoadConfig::getParam( $rH_cfg, "default", 'sampleOutputRoot' ) . $sampleName . '/output_jobs/' . "\n";
}

sub printSubmitCmd {
    my $rH_cfg         = shift;
    my $stepName       = shift;
    my $jobNameSuffix  = shift;
    my $jobIdPrefix    = shift;
    my $dependancyName = shift;
<<<<<<< HEAD
    my $sampleName = shift;
    my $command = shift;
    my $workDirectory = shift;
=======
    my $sampleName     = shift;
    my $command        = shift;
>>>>>>> 26201201869edf380313f3e0dbe689139bb0019a

    my $jobIdVarName = uc($jobIdPrefix) . '_JOB_ID';

<<<<<<< HEAD
    if(!(defined $workDirectory)){
      $workDirectory = '`pwd`';
      if(LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "true") {
        $command =~ s/\\\$/\\\\\\\$/g;
        print $jobIdVarName.'=`';
        $workDirectory = '\`pwd\`';
      }
    }
    print 'echo "'.$command.'" | ';
    print LoadConfig::getParam($rH_cfg, $stepName, 'clusterSubmitCmd');
    print ' ' . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOtherArg');
    print ' ' . LoadConfig::getParam($rH_cfg, $stepName, 'clusterWorkDirArg') . $workDirectory;
    print ' ' . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOutputDirArg') . $sampleName.'/output_jobs/';
    my $jobName = $stepName.'.'.$sampleName;
    if(defined($jobNameSuffix) && length($jobNameSuffix) > 0) {
      $jobName .= '.'.$jobNameSuffix;
=======
    my $currentDirCommand = '`pwd`';
    if ( LoadConfig::getParam( $rH_cfg, $stepName, 'clusterCmdProducesJobId' ) eq "true" ) {
        $command =~ s/\\\$/\\\\\\\$/g;
        print $jobIdVarName. '=`';
        $currentDirCommand = '\`pwd\`';
    }
    print 'echo "' . $command . '" | ';
    print LoadConfig::getParam( $rH_cfg, $stepName, 'clusterSubmitCmd' );
    print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterOtherArg' );
    print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterWorkDirArg' ) . ' ' . $currentDirCommand;
    print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterOutputDirArg' ) . ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'sampleOutputRoot' ) . $sampleName . '/output_jobs/';
    my $jobName = $stepName . '.' . $sampleName;
    if ( defined($jobNameSuffix) && length($jobNameSuffix) > 0 ) {
        $jobName .= '.' . $jobNameSuffix;
>>>>>>> 26201201869edf380313f3e0dbe689139bb0019a
    }
    print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterJobNameArg' ) . ' ' . $jobName;
    print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterWalltime' );
    print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterQueue' );
    print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterCPU' );
    if ( defined($dependancyName) ) {
        print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterDependencyArg' ) . $dependancyName;
    }
    print ' ' . LoadConfig::getParam( $rH_cfg, $stepName, 'clusterSubmitCmdSuffix' );
    if ( LoadConfig::getParam( $rH_cfg, $stepName, 'clusterCmdProducesJobId' ) eq "true" ) {
        print '`';
    }
    print "\n\n";

    if ( LoadConfig::getParam( $rH_cfg, $stepName, 'clusterCmdProducesJobId' ) eq "false" ) {
        print $jobIdVarName. '=' . $jobName . "\n";
    }

    return $jobIdVarName;
}
1;
