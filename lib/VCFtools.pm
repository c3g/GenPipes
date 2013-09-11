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

package VCFtools;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub annotateMappability {
    my $rH_cfg      = shift;
    my $inputVCF    = shift;
    my $outputVCF   = shift;

    my $outDate = -M $outputVCF;
    my $inDate = -M $inputVCF;
  
    my $command;

    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'annotateMappability', 'moduleVersion.vcftools').' ;';
        $command .= ' module load '.LoadConfig::getParam($rH_cfg, 'annotateMappability', 'moduleVersion.tabix').' ;';
        $command .= ' vcf-annotate -d \"key=INFO,ID=MIL,Number=1,Type=String,Description='."'".'Mappability annotation. 300IS 40SD 1SHI. HC = to high coverage (>400), LC = to high coverage (<50), MQ = to low mean mapQ (<20), ND = no data at the position'."'".'\"';
        $command .= ' -c CHROM,FROM,TO,INFO/MIL ';
        $command .= ' -a '.LoadConfig::getParam($rH_cfg, 'annotateMappability', 'referenceMappabilityBedIndexed');
        $command .= ' '.$inputVCF;
        $command .= ' > '.$outputVCF;
    #}
    return $command;
}

sub indexVCF {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $inputVCF    = shift;
    my $outputVCF   = shift;

    my $outDate = -M $outputVCF;
    my $inDate = -M $inputVCF;
  
    my $command;

    # -M gives modified date relative to now. The bigger the older.
    #if(!defined($outDate) || !defined($inDate) || $inDate < $outDate) {
        $command .= 'module load '.LoadConfig::getParam($rH_cfg, 'indexVCF', 'moduleVersion.tabix').' ;';
        $command .= ' bgzip -c';
        $command .= ' '.$inputVCF;
        $command .= ' > '.$outputVCF;
        $command .= ' ; ';
        $command .= ' tabix -p vcf -f';
        $command .= ' '.$outputVCF;
    #}
    return $command;
}
1;
