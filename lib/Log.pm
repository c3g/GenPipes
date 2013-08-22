#!/usr/env/perl

=head1 NAME

I<Logh>

=head1 SYNOPSIS

Log-> readRestartFile()

=head1 DESCRIPTION

B<Log> is a library to manage the accessor of internal log file and to prodiu



=head1 AUTHOR
B<Joel Fillon> - I<TO_BE_ADDED>
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package Log;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------

sub readRestartFile {
  my $rH_cfg      = shift;
  my $restartLogFilePath  = shift;
  my $rH_restartList = shift;
  
  
  
}

1;
