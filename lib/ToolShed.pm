#!/usr/env/perl

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

package ToolShed;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Cwd 'abs_path';
use File::Basename;

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

    my $outDate = -M $outputVCF;
    my $inDate = -M $inputVCF;
  
    my $toolShedDir = getToolShedDir();

    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= $toolShedDir.'/filterLongIndel.pl ';
        $command .= ' '.$inputVCF;
        $command .= ' > '.$outputVCF;
    #}
    return $command;
}

1;
