#!/usr/bin/env perl

=head1 NAME

I<SVTools>

=head1 SYNOPSIS

Picard->merge()

=head1 DESCRIPTION

B<SVTools> is a library to analyse BAMs for Structural Variants

=head1 AUTHOR

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package Tools;

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
use Cwd 'abs_path';
use File::Basename;
use LoadConfig;
use Socket;
use Sys::Hostname;

# SUB
#-----------------------
sub getToolShedDir {
  my $currentDir = dirname(__FILE__);
  return abs_path($currentDir.'/../tool_shed');
}

sub filterNStretches {
  my $rH_cfg      = shift;
  my $sampleName  = shift;
  my $inputVCF    = shift;
  my $outputVCF   = shift;

  my $toolShedDir = getToolShedDir();

  my $rO_job = new Job();
  $rO_job->testInputOutputs([$inputVCF], [$outputVCF]);

  if (!$rO_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['metrics' , 'moduleVersion.tools']]) . ' &&';
    $command .= ' \$PERL_TOOLS/filterLongIndel.pl ';
    $command .= ' ' . $inputVCF;
    $command .= ' > ' . $outputVCF;

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub generateIntervalList {
  my $rH_cfg   = shift;
  my $dict     = shift;
  my $bedFile  = shift;
  my $output   = shift;

  my $rO_job = new Job();
  if(!defined($dict)) {
    $dict = LoadConfig::getParam($rH_cfg, 'default', 'referenceSequenceDictionary', 1, 'filepath');
  }

  $rO_job->testInputOutputs([$dict, $bedFile], [$output]);

  if (!$rO_job->isUp2Date()) {
    my $command;
    $command .= LoadConfig::moduleLoad($rH_cfg, [['default' , 'moduleVersion.tools'], ['default' , 'moduleVersion.perl']]) . ' &&';
    $command .= ' bed2IntervalList.pl';
    $command .= ' --dict ' . $dict;
    $command .= ' --bed ' . $bedFile;
    $command .= ' > ' . $output;

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub mugqicLog {
  my $pipeline = shift;
  my $steps = shift;
  my $samples = shift;

  my $server = "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi";
  my $hostname = hostname;
  # Retrieve client local IP address
  my $ip = inet_ntoa((gethostbyname(hostname))[4]);

  my $request =
    "hostname=$hostname&" .
    "ip=$ip&" .
    "pipeline=$pipeline&" .
    "steps=$steps&" .
    "samples=$samples"
  ;

  print "#" . "-" x 79 . "\n";
  print "# Call home with pipeline statistics\n";
  print "#" . "-" x 79 . "\n";
  print "wget \"$server?$request\" --quiet --output-document=/dev/null\n";
}

1;
