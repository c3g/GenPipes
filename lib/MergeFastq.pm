#!/usr/env/perl

=head1 NAME

I<MergeFastq> is a lib that merges calls a python
script which merges fastq files if needed.

=head1 SYNOPSIS

MergeFastq::mergeFiles(%ref_hash_config, $runType, $fileAlias, $outDir)

$fileAlias is a file containing the fastaq information (path, genome name ...)

$runType: 1 pair end

           0 single end

$outDir is the dirname of FileAlias.
		   
=head1 DESCRIPTION

B<MergeFastq> is a library Merges fastq files before assembly


=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package MergeFastq;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Data::Dumper;
use Config::Simple;
use File::Basename;
use Cwd 'abs_path';

sub mergeFiles {
    my $rH_cfg    = shift;
    my $runType   = shift;
    my $fileAlias = shift;
    my $outDir    = dirname($fileAlias);

    my $command = '';
    my %retVal;

    $runType = ( $runType eq "PAIRED_END" ) ? 1 : 0;

    my $mergeStatus = `$rH_cfg->{'merge.python'} -f $fileAlias -d $runType -o $outDir`;
    if ( $mergeStatus == 0 ) {

        #---------------------------------------
        # Todo
        # Need to Implement command line merge
        #--------------------------------------
        $command = ' pyhton Do a merge';

    }
    else {
        print "No Merge step necessary\n";
        $command = 'No merge';

    }
    $retVal{'command'} = $command;
    return ( \%retVal );

}

1;

