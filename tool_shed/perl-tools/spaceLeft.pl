#!/usr/bin/perl  

use Filesys::Df;
use Getopt::Long;
use Net::SMTP;
use strict;

&main();

sub main {
  my @testDirectories;
  my @to;
  my $verbose;
  my $fromEmail;
  my $result = GetOptions ( "from=s" => \$fromEmail,
                            "to=s" => \@to,
                            "dir=s" => \@testDirectories,
                            "verbose" => \$verbose);

	my $message = "";
	for my $directory (@testDirectories) {
		my $dirInfo = df($directory, 1);
    if(defined($directory)){
      my $dirSize = $dirInfo->{"bfree"};
			if($dirSize < 4 * 1024 * 1024 * 1024 * 1024) {  #4TB
				for (qw(B KB MB GB TB PB)) {
					if($dirSize < 1024) {
        	  $message .= sprintf("%s Space left: %.1f%s\n", $directory, $dirSize, $_);
						last;
					}
					$dirSize=$dirSize/1024;
				}
			}
		}
		else {
			warn "Undefined $directory\n";
		}
	}

	if(length($message) > 0) {
		sendEmail($fromEmail,\@to, $message);
    print $message;
	}
}

sub sendEmail {
  my $fromEmail = shift;
  my $rA_to = shift;
  my $message = shift;

  my $smtp = Net::SMTP->new('mailhost.mcgill.ca');
  if(!defined($smtp) || !($smtp)) {
    print STDERR "SMTP ERROR: Unable to open smtp session.\n";
    return 1;
  }

  if (! ($smtp->mail($fromEmail) ) ) {
    print STDERR "SMTP ERROR: Cannot set from: ".$fromEmail.".\n";
    return 2;
  }

  if (! ($smtp->to(@$rA_to) ) ) {
    print STDERR "SMTP ERROR: Cannot set to: ".join(',', @$rA_to).".\n";
    return 3;
  }

  my $tos = join(',',@$rA_to);

  my $subject = "Subject: Make space for HiSeqs\n";

  $smtp->data();
  $smtp->datasend("To: ".$tos."\n");
  $smtp->datasend($subject);
  $smtp->datasend("\n");
	$smtp->datasend($message);
  $smtp->dataend();

  $smtp->quit;

  return 0;
}

