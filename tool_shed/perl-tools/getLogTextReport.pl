#!/usr/bin/perl -w

# Strict Pragmas
#---------------------
use strict;
use warnings;
#---------------------

use File::Basename;
use Log;

# Check and assign the parameters
$ARGV[0] or die "Usage: perl ".basename($0)." <job_log_list_file>\n";
my $jobLogList = $ARGV[0];

print Log::getLogTextReport($jobLogList);
