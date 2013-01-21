#!/usr/env/perl

=head1 NAME

I<LoadModules>

=head1 SYNOPSIS

LoadModules->getModules(file_name)

=head1 DESCRIPTION

B<LoadModules> is a library that reads a list of modules (unix modulecmd)
and returns a 'module add PACKAGE' command. The command can be used directly on
a PBS file or looped and add to a bash script.

Input = file_name

Output = array


=head1 AUTHOR

B<David Morais> - I<dmorais@cs.bris.ac.uk>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package  LoadModules;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Data::Dumper;
use Pod::Usage;

# SUB
#-----------------------
sub getModules {
    my ( $self, $file ) = @_;
    my @modules;

    open( IN, $file ) || die "could not open $file, $!\n";
    while (<IN>) {
        chomp;
        push( @modules, "module add $_\n" );

    }
    return @modules;
}

1;
