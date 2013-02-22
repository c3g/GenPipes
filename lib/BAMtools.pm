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

package BAMtools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub countBins {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $tumorBam    = shift;
    my $window      = shift;
    my $normType    = shift;
    my $outputFile  = shift;
    my $normalBam   = shift; # can be undef for non paired run

    my $outDate = -M $outputFile;
    my $inDate = -M $tumorBam;
  
    my $command;
    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'countBins', 'moduleVersion.bamtools').' ; ';
        $command .= ' java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'countBins', 'tmpDir').' '.LoadConfig::getParam($rH_cfg, 'countBins', 'extraJavaFlags');
        $command .= ' -Xmx1500M -jar \${BAMTOOLS_JAR} bincounter';
        $command .= ' --norm '.$normType;
        $command .= ' --minMapQ '.LoadConfig::getParam($rH_cfg, 'countBins', 'minMapQ');
        $command .= ' --bam '.$tumorBam;
        if(defined($normalBam)) {
            $command .= ' --refbam '.$normalBam;
        }
        $command .= ' --window '.$window;
        $command .= ' > '.$outputFile;
    #}
    return $command;
}
1;
