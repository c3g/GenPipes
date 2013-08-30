#!/usr/env/perl

=head1 NAME

I<DiffExpression>

=head1 SYNOPSIS

DiffExpression:sub(args)

B<DiffExpression:stats>(%ref_hash_config, $sample_name, %ref_hash_laneInfo, $group) 

All subroutines return a ref_hash with the command line


=head1 DESCRIPTION

B<DiffExpression> is a library to create differencial expression analysis.

=head1 AUTHOR

David Morais dmorais@cs.bris.ac.uk
B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug


=cut

package DiffExpression;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use PipelineUtils;
use Data::Dumper;
use Config::Simple;

our $rH_cfg;
our $sampleName;
our $rH_laneInfo;
our $readFile;
our $group;

sub edger {
    $rH_cfg      = shift;
    $sampleName  = shift;
    $rH_laneInfo = shift;
    $group       = shift;
    my %retVal;

    my $laneDirectory = 'DGE/' . $group . "/";
    my $up2date = PipelineUtils::testInputOutputs([$laneDirectory . 'matrix.csv '], [$laneDirectory]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= ' module add gcc/4.7.0 ; module add R ;';
        $command .= ' Rscript ' . $rH_cfg->{'diffExpress.edger'} . ' -d ' . $rH_cfg->{'diffExpress.designFile'};
        $command .= ' -c ' . $laneDirectory . 'matrix.csv ';
        $command .= ' -o ' . $laneDirectory . ' ;';
        $command .= ' Rscript ' . $rH_cfg->{'diffExpress.deseq'} . ' -d ' . $rH_cfg->{'diffExpress.designFile'};
        $command .= ' -c ' . $laneDirectory . 'matrix.csv ';
        $command .= ' -o ' . $laneDirectory . ' ;';
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }

    return $ro_job;

}

sub edgerPortable {
      my $rH_cfg        = shift;
      my $designFile    = shift;
      my $countMatrix   = shift;
      my $outputDir     = shift;

    my $up2date = PipelineUtils::testInputOutputs([$countMatrix], [$outputDir]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'diffExpress','moduleVersion.tools') .' ' .LoadConfig::getParam($rH_cfg, 'diffExpress','moduleVersion.cranR') .' ;';
        $command .= ' Rscript \$R_TOOLS/edger.R -d ' .$designFile;
        $command .= ' -c ' .$countMatrix;
        $command .= ' -o ' .$outputDir;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }

    return $ro_job;
}

sub deseq {
    my $rH_cfg        = shift;
    my $designFile    = shift;
    my $countMatrix   = shift;
    my $outputDir     = shift;
    
    my $up2date = PipelineUtils::testInputOutputs([$countMatrix], [$outputDir]);
    my $ro_job = new Job(!defined($up2date));

    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'diffExpress','moduleVersion.tools') .' ' .LoadConfig::getParam($rH_cfg, 'diffExpress','moduleVersion.cranR') .' ;';
        $command .= ' Rscript \$R_TOOLS/deseq.R -d ' .$designFile;
        $command .= ' -c ' .$countMatrix;
        $command .= ' -o ' .$outputDir;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }

    return $ro_job;
}

sub goseq {
    my $rH_cfg        = shift;
    my $resultFile    = shift;
    my $outputFile    = shift;
    my $columns       = shift;
    my $method        = shift;
    
    my $maxResult = LoadConfig::getParam($rH_cfg, 'goseq','maxGoResult');
    my $geneSizeFile = LoadConfig::getParam($rH_cfg, 'goseq','geneSizeFile');
    my $goLinkFile = LoadConfig::getParam($rH_cfg, 'goseq','goLinkFile');
    my $geneIdType = LoadConfig::getParam($rH_cfg, 'goseq','geneIdType');
    my $option = '';
    if (defined($maxResult) && $maxResult ne "" && $maxResult ne "0") {
        $option = ' -m ' .$maxResult;
    }
    if (defined($geneSizeFile) && $geneSizeFile ne "") {
        $option .= ' -a ' .$geneSizeFile;
    }
    if (defined($goLinkFile) && $goLinkFile ne "") {
        $option .= ' -G ' .$goLinkFile;
    }
    if (defined($geneIdType) && $geneIdType ne "") {
        $option .= ' -i ' .$geneIdType;
    }
    
    my $up2date = PipelineUtils::testInputOutputs([$resultFile], [$outputFile]);
    my $ro_job = new Job(!defined($up2date));
    
    if (!$ro_job->isUp2Date()) {
        my $command;
        $command .= 'module load ' .LoadConfig::getParam($rH_cfg, 'goseq','moduleVersion.tools') .' ' .LoadConfig::getParam($rH_cfg, 'goseq','moduleVersion.cranR') .' &&';
        $command .= ' Rscript \$R_TOOLS/goseq.R -d ' .$resultFile;
        $command .= ' -c ' .$columns;
        $command .= ' -t ' .LoadConfig::getParam($rH_cfg, 'goseq','goAnnotation');
        $command .= ' -k ' .LoadConfig::getParam($rH_cfg, 'goseq','referenceEnsemble2symbol');
        $command .= ' -s ' .LoadConfig::getParam($rH_cfg, 'goseq','referenceUCSCname');
        $command .= $option;
        $command .= ' -o ' .$outputFile;
        $command .= ' ' . $up2date;

        $ro_job->addCommand($command);
    }

    return $ro_job;
}

1;

