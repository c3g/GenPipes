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
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'countBins', 'moduleVersion.java').' ' .LoadConfig::getParam($rH_cfg, 'countBins', 'moduleVersion.bamtools').' ; ';
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

sub deleteDuplicates {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $pair1           = shift;
    my $pair2           = shift;
    my $single          = shift;
    my $optOutputPrefix = shift;


    my $rH_retVal;
    if (defined($pair1) && defined($pair2)) {
        $rH_retVal = deletePairedDuplicates($rH_cfg, $sampleName, $pair1, $pair2, $optOutputPrefix);
    }
    elsif (defined($single)) {
        $rH_retVal = deleteSingleDuplicates($rH_cfg, $sampleName, $single, $optOutputPrefix);
    }
    else {
        die "Unknown runType. \n";
    }

    return $rH_retVal;
}

sub deletePairedDuplicates {
    my $rH_cfg = shift;
    my $sampleName = shift;
    my $pair1           = shift;
    my $pair2           = shift;
    my $outputPrefix    = shift;
    my %retVal;
	
    my $command             = '';
    my $outputFastqPair1Name = $outputPrefix . '.pair1.dup.fastq.gz';
    my $outputFastqPair2Name = $outputPrefix . '.pair2.dup.fastq.gz';
    my $pair1FileDate = -M $pair1;
    my $pair2FileDate = -M $pair2;

    my $currentFileDate = $pair2FileDate;
    if ( defined($pair1FileDate) && $pair1FileDate > $pair2FileDate ) {
        $currentFileDate = $pair1FileDate;
    }
    my $outDate = -M $outputFastqPair1Name;
	
    if ( !defined($currentFileDate) || !defined($outDate) || $currentFileDate > -M $outDate ){
        $command .= 'module add '. LoadConfig::getParam($rH_cfg, 'default','moduleVersion.java') ;
        $command .= ' ' . LoadConfig::getParam($rH_cfg, 'duplicate','moduleVersion.bamtools') . ';' ;
        $command .= ' java '.LoadConfig::getParam($rH_cfg, 'duplicate', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'duplicate', 'dupRam').' -jar  \$BAMTOOLS_JAR filterdups' .  ' --read1 ' . $pair1 . ' --read2 ' . $pair2;
        $command .= ' -k 20 -o 15;';
        $command .= ' mv ' . $pair1 . '.dup.read1.gz ' . $outputFastqPair1Name . ';';
        $command .= ' mv ' . $pair2 . '.dup.read2.gz ' . $outputFastqPair2Name . ';';
    }
    $retVal{'command'} = $command;
    $retVal{'pair1'}   = $outputFastqPair1Name;
    $retVal{'pair2'}   = $outputFastqPair2Name;

    return ( \%retVal );

}

sub deleteSingleDuplicates {
    my $rH_cfg = shift;
    my $sampleName = shift;
    my $single          = shift;
    my $outputPrefix    = shift;
    my %retVal;
	
    my $command         = '';
    my $outputFastqName = $outputPrefix . '.single.dup.fastq.gz';
    my $currentFileDate = -M $outputFastqName;

    my $outDate = -M $outputFastqName;
    if ( !defined($currentFileDate) || !defined($outDate) || $currentFileDate > -M $outputFastqName) {
        $command .= 'module add '. LoadConfig::getParam($rH_cfg, 'default','moduleVersion.java') ;
        $command .= ' ' . LoadConfig::getParam($rH_cfg, 'duplicate','moduleVersion.bamtools') . ';' ;
        $command .= ' java '.LoadConfig::getParam($rH_cfg, 'duplicate', 'extraJavaFlags').' -Xmx'.LoadConfig::getParam($rH_cfg, 'duplicate', 'dupRam').' -jar \$BAMTOOLS_JAR filterdups' . ' --read1 ' . $single;
        $command .= ' -k 20 -o 15;';
        $command .= ' mv ' . $single . '.dup.read1.gz ' . $outputFastqName . ';';
    }
    $retVal{'command'} = $command;
    $retVal{'single1'} = $outputFastqName;

    return ( \%retVal );
}

1;

