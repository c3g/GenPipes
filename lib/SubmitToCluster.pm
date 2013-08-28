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
use Scalar::Util 'blessed';
use Cwd;
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
    my $sampleName     = shift;
    my $command        = shift;
    my $outputDir      = shift;
    my $workDirectory  = shift;
    my $commandIdx     = shift;

    my $isBlessed = defined(blessed( $command ));
    my $ro_job = undef;
    if($isBlessed) {
        $ro_job = $command;

        if($ro_job->isUp2Date()) {
          return undef;
        }

        if(!defined($commandIdx)) {
          $commandIdx = 0;
        }
        $command = $ro_job->getCommand($commandIdx);
    }

    my $jobIdVarName = uc( $jobIdPrefix ) . '_JOB_ID';
    $jobIdVarName =~ s/\W/_/g;
    #$jobIdVarName = ~ s/^[A-Za-z0-9\_]/_/g;

    ### TO DO modify the output dir to be more portable

    if(!(defined $workDirectory)){
      $workDirectory = getcwd();
    }

    if(!(defined $outputDir)){
      $outputDir = $sampleName;
    }
    if(LoadConfig::getParam($rH_cfg, $stepName, 'clusterCmdProducesJobId') eq "true") {
        print $jobIdVarName.'=$(';
    }
    print 'echo "'.$command.'" | ';
    print LoadConfig::getParam($rH_cfg, $stepName, 'clusterSubmitCmd');
    print ' ' . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOtherArg');
    print ' ' . LoadConfig::getParam($rH_cfg, $stepName, 'clusterWorkDirArg') . ' ' . $workDirectory;
    print ' ' . LoadConfig::getParam($rH_cfg, $stepName, 'clusterOutputDirArg') .' '  .$outputDir .'/output_jobs/';
    my $jobName ;
    if(defined($sampleName) && length($sampleName) > 0) {
        $jobName = $stepName.'.'.$sampleName;
    }
    else {
        $jobName = $stepName;
    }
    if(defined($jobNameSuffix) && length($jobNameSuffix) > 0) {
      $jobName .= '.'.$jobNameSuffix;
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
        print ')';
    }
    print "\n\n";

    if ( LoadConfig::getParam( $rH_cfg, $stepName, 'clusterCmdProducesJobId' ) eq "false" ) {
        print $jobIdVarName. '=' . $jobName . "\n";
    }

    if($isBlessed) {
      $ro_job->setCommandJobId($commandIdx, '$'.$jobIdVarName);
    }
    return $jobIdVarName;
}
1;
