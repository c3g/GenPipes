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
use LWP::UserAgent;
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

sub mugqicLog {
  my $pipeline = shift;
  my $steps = shift;
  my $samples = shift;

  my $userAgent = LWP::UserAgent->new;
  my $server = "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi";
  my $request = HTTP::Request->new(POST => $server);

  my $hostname = hostname;
  # Retrieve client local IP address
  my $ip = inet_ntoa((gethostbyname(hostname))[4]);

  $request->content(
    "hostname=$hostname&" .
    "ip=$ip&" .
    "pipeline=$pipeline&" .
    "steps=$steps&" .
    "samples=$samples"
  );

  my $response = $userAgent->request($request);
  if ($response->is_success) {
    my $message = $response->decoded_content;
    print STDERR "MUGQIC remote log sent successfully. $message\n";
  } else {
    print STDERR "MUGQIC remote log failed to be sent.\n";
    print STDERR "HTTP GET error code: ", $response->code, "\n";
    print STDERR "HTTP GET error message: ", $response->message, "\n";
  }
}

1;
