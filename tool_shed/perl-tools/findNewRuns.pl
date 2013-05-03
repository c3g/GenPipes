#!/usr/bin/perl

require 5.008;

use strict;
use POSIX qw( strftime );
use Net::SMTP;
use XML::Simple;
use Data::Dumper;
#use Time::localtime;
#use File::stat;
use Getopt::Long;

my $version = "0.1";

&main();

sub main {
  my @runDirectories;
  my @to;
  my $verbose;
  my $fromEmail;
  my $result = GetOptions ( "from=s" => \$fromEmail,
                            "to=s" => \@to,
                            "runDir=s" => \@runDirectories,
                            "verbose" => \$verbose);

  my $nanuqRunDoneFileName = "nanuqRunFinished-seen.empty";
  my @newRunsDone;
  my @newRunsDoneFiles;
  my @oldRunsDone;
  for my $rootDir (@runDirectories) {
    if(! -d $rootDir) {
      warn "'".$rootDir."' is not a directory. Ignored\n";
    }

    opendir(ROOT_DIR, $rootDir) or die "Couldn't open directory ".$rootDir."\n";
    my @rootFiles =  grep { /^[^\.]/ && -d "$rootDir/$_" } readdir(ROOT_DIR);
    closedir(ROOT_DIR);

    for my $partialRunDir (@rootFiles) {
      my $runDir = $rootDir.'/'.$partialRunDir;
      if(! -d $runDir) {
        next;
      }

			if(! -e $runDir."/RunInfo.xml") {
				next;
			}

			my $xml = new XML::Simple ('ForceArray' => [ 'Read' ]);
			my $runInfo = $xml->XMLin($runDir."/RunInfo.xml");
			my $nbReads = @{$runInfo->{'Run'}->{'Reads'}->{'Read'}};
			my $fileToLookFor;
			my $isIndexed= "N";
			if($nbReads > 1) {
				$isIndexed = $runInfo->{'Run'}->{'Reads'}->{'Read'}->[1]->{'IsIndexedRead'};
			}
			else {
				$isIndexed = $runInfo->{'Run'}->{'Reads'}->{'Read'}->[0]->{'IsIndexedRead'};
			}
			$fileToLookFor="Basecalling_Netcopy_complete_Read".$nbReads.".txt";

      opendir(RUN_DIR, $runDir) or die "Couldn't open possible run directory ".$runDir."\n";
      my @runFiles =  grep { /^[^\.]/ && -f "$runDir/$_" } readdir(RUN_DIR);
      closedir(RUN_DIR);

      my $runDone = 0;
      my $oldRun = 0;
      for my $partialRunFile (@runFiles) {
        #my $runFile = $runDir.'/'.$partialRunFile;
        if($partialRunFile eq $fileToLookFor) {
          $runDone = 1;
        } elsif ($partialRunFile eq $nanuqRunDoneFileName) {
          $oldRun = 1;
        }
      }

      if($runDone) {
        $runDir =~ /.*\/(\d{6}_[^_]+_[0-9]{4}_[^\/]+(_.*)?)/;
        my $completeFile = $runDir.'/'.$fileToLookFor;
        my $formattedCreationTime = POSIX::strftime( "%Y-%m-%d %H:%M:%S", localtime( (stat($completeFile))[9] ) );
        if($oldRun) {
					my $message = "Old: ".$runDir." -- IsIndexed: ".$isIndexed." -- ".$formattedCreationTime."\n";
          if($verbose) {
						print STDERR $message;
					}
          push(@oldRunsDone, $message);
        }
        else {
					my $message = "New: ".$runDir." -- IsIndexed: ".$isIndexed." -- ".$formattedCreationTime."\n";
          if($verbose) {
						print STDERR $message;
					}
          push(@newRunsDone, $message);
          push(@newRunsDoneFiles, $runDir);
          
          my $command = '/sb/programs/analyste/software/tools/processRuns.pl --nanuqAuthFile $HOME/.nanuqAuth.txt --email '.$fromEmail.' --runDir '.$runDir.' > '.$runDir.'/run.sh && chmod ug+rwx '.$runDir.'/run.sh';
          print STDERR "Command: $command\n";
          system($command);
        }
      }
    }
  }

  my $errorCode = 0;
  if(@newRunsDone > 0) {
    $errorCode = sendEmail($fromEmail, \@to, \@newRunsDone, \@oldRunsDone);
    if($errorCode == 0) {
      for my $newRun (@newRunsDoneFiles) {
        my $touchFile = $newRun."/".$nanuqRunDoneFileName;
        open(TOUCH, ">", $touchFile) or warn "Cannot touch: ".$touchFile."\n";
				close(TOUCH);
        if($verbose) {
					print STDERR "Touched: ".$touchFile."\n";
				}
      }
    }
  }
  exit $errorCode;
}

sub sendEmail {
  my $fromEmail = shift;
  my $rA_to = shift;
  my $rA_newRunsDone = shift;
  my $rA_oldRunsDone = shift;

  my $smtp = Net::SMTP->new('192.168.32.1');
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

	my $subject = "Subject: Runs done. New: ";
  $subject .= @$rA_newRunsDone;
  $subject .= " Old: ";
  $subject .= @$rA_oldRunsDone;
  $subject .= "\n";

  $smtp->data();
  $smtp->datasend("To: ".$tos."\n");
  $smtp->datasend($subject);
  $smtp->datasend("\n");
  for my $newRun (@$rA_newRunsDone) {
    $smtp->datasend($newRun);
  }
  $smtp->datasend("\n");
  for my $oldRun (@$rA_oldRunsDone) {
    $smtp->datasend($oldRun);
  }
  $smtp->dataend();

  $smtp->quit;

  return 0;
}

