#!/usr/env/perl

=head1 NAME

I<SmrtAnalysis>

=head1 SYNOPSIS

SmrtAnalysis->run()

=head1 DESCRIPTION

B<Picard> This a library to analyze PacBio data using the SmrtAnalysis suite.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package SmrtAnalysis;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin";

# Dependencies
#-----------------------
use Job;
use LoadConfig;

# SUB
#-----------------------
sub run {
  my $rH_cfg     = shift;
  my $paramsXml  = shift;
  my $inputXml   = shift;
  my $outdir     = shift;
  my $log        = shift;
  my $tmpdir     = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$paramsXml, $inputXml], [$log]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' smrtpipe.py';
    $cmd .= ' -D NPROC=' . LoadConfig::getParam($rH_cfg, 'smrtanalysis', 'num_threads');
    $cmd .= ' -D TMP=' . $tmpdir;
    $cmd .= ' --params=' . $paramsXml;
    $cmd .= ' --output=' . $outdir;
    $cmd .= ' --debug'; # ... Leave it there for now.
    $cmd .= ' xml:' . $inputXml;
    $cmd .= ' > ' . $log . ' 2>&1';

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub filtering {
  my $rH_cfg         = shift;
  # Fofn
  my $fofn           = shift;
  my $infilesXml     = shift;
  my $outFofn        = shift;
  # xml
  my $refParamsXml   = shift;
  my $currParamsXml  = shift;
  # smrtpipe
  my $outdir         = shift;
  my $log            = shift;
  my $tmpdir         = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$fofn], [$outdir . '/data/filtered_subreads.fastq']);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    # Fofn
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['prinseq', 'moduleVersion.prinseq'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' fofnToSmrtpipeInput.py ';
    $cmd .= $fofn;
    $cmd .= ' > ' . $infilesXml;
    $cmd .= ' && cp ' . $fofn . ' ' . $outFofn;
    $cmd .= ' && ';
    # Xml
    $cmd .= 'memtime';
    $cmd .= " sed -e \'s/MINSUBREADLENGTH/" . LoadConfig::getParam($rH_cfg, 'filtering', 'minSubReadLength') . "/g\'";
    $cmd .= " -e \'s/MINREADLENGTH/" . LoadConfig::getParam($rH_cfg, 'filtering', 'minReadLength') . "/g\'";
    $cmd .= " -e \'s/MINQUAL/" . LoadConfig::getParam($rH_cfg, 'filtering', 'minQual') . "/g\' < " . $refParamsXml . " > " . $currParamsXml;
    $cmd .= ' && ';
    # SmrtPipe
    $cmd .= 'memtime';
    $cmd .= ' smrtpipe.py';
    $cmd .= ' -D NPROC=' . LoadConfig::getParam($rH_cfg, 'filtering', 'num_threads');
    $cmd .= ' -D TMP=' . $tmpdir;
    $cmd .= ' --params=' . $currParamsXml;
    $cmd .= ' --output=' . $outdir;
    $cmd .= ' --debug';
    $cmd .= ' xml:' . $infilesXml;
    $cmd .= ' > ' . $log;# . ' 2>&1';
    # Fasta to Fastq
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' prinseq-lite.pl';
    $cmd .= ' -verbose';
    $cmd .= ' -fastq ' . $outdir . '/data/filtered_subreads.fastq';
    $cmd .= ' -out_format 1';
    $cmd .= ' -out_good ' . $outdir . '/data/filtered_subreads';
    #$cmd .= ' > /dev/null 2>&1';

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}


sub fofns {
  my $rH_cfg     = shift;
  my $fofn       = shift;
  my $outputXml  = shift;
  my $outFofn    = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$fofn], [$outputXml]);

  if (!$ro_job->isUp2Date()) {
  my $cmd = '';

    # Then do smrtpipe to generate input.xml
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' fofnToSmrtpipeInput.py ';
    $cmd .= $fofn;
    $cmd .= ' > ' . $outputXml;
    $cmd .= ' && cp ' . $fofn . ' ' . $outFofn;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub referenceUploader {

  my $rH_cfg     = shift;
  my $prefix     = shift;
  my $sampleName = shift;
  my $fasta      = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$fasta], [$prefix . "/reference.info.xml"]);

  if (!$ro_job->isUp2Date()) {
  my $cmd = '';

    # Preload assembled contigs as reference
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh && ';
    $cmd .= ' memtime';
    $cmd .= ' referenceUploader';
    $cmd .= ' -c';
    $cmd .= ' -p ' . $prefix;
    $cmd .= ' -n ' . $sampleName;
    $cmd .= ' -f ' . $fasta;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;

}

sub runCommand {

  my $rH_cfg     = shift;
  my $string     = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs(undef, undef );

  if (!$ro_job->isUp2Date()) {
  my $cmd;

    # Preload assembled contigs as reference
    $cmd = $string;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;

}

sub blasr {
  my $rH_cfg      = shift;
  my $infile      = shift;
  my $infileLong  = shift;
  my $outfile     = shift;
  my $outfileFofn = shift;
  my $sam         = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile, $infileLong], [$outfile, $outfileFofn]);

  if (!$ro_job->isUp2Date()) {
    #blasr filtered_subreads.fa filtered_longreads -out seeds.m4 -m 4 -nproc 8 -bestn 24 -nCandidates 24 -noSplitSubreads -minReadLength 200 -maxScore -1000 -maxLCPLength 16
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['default', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' blasr';
    $cmd .= ' ' . $infile;
    $cmd .= ' ' . $infileLong;
    $cmd .= ' -out ' . $outfile;
    $cmd .= ' -m ' . LoadConfig::getParam($rH_cfg, 'blasr', 'm');
    $cmd .= ' -nproc ' . LoadConfig::getParam($rH_cfg, 'blasr', 'num_threads');
    $cmd .= ' -bestn ' . LoadConfig::getParam($rH_cfg, 'blasr', 'bestn');
    $cmd .= ' -nCandidates ' . LoadConfig::getParam($rH_cfg, 'blasr', 'nCandidates');
    $cmd .= ' -noSplitSubreads';
    $cmd .= ' -minReadLength ' . LoadConfig::getParam($rH_cfg, 'blasr', 'minReadLength');
    $cmd .= ' -maxScore ' . LoadConfig::getParam($rH_cfg, 'blasr', 'maxScore');
    $cmd .= ' -maxLCPLength ' . LoadConfig::getParam($rH_cfg, 'blasr', 'maxLCPLength');
    if ($sam eq "sam") {
      $cmd .= ' -sam';
    }
    $cmd .= ' && echo ' . $outfile . ' > ' . $outfileFofn;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub m4topre {
  my $rH_cfg     = shift;
  my $infile     = shift;
  my $allm4      = shift;
  my $subreads   = shift;
  my $outfile    = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfile]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' m4topre.py';
    $cmd .= ' ' . $infile;
    $cmd .= ' ' . $allm4;
    $cmd .= ' ' . $subreads;
    $cmd .= ' ' . LoadConfig::getParam($rH_cfg, 'm4topre', 'bestn');
    $cmd .= ' > ' . $outfile;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub pbdagcon {
  my $rH_cfg       = shift;
  my $infile       = shift;
  my $outfile      = shift;
  my $outfileFastq = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfile, $outfileFastq]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh && ';
    $cmd .= ' memtime';
    $cmd .= ' pbdagcon';
    $cmd .= ' -a ';
    $cmd .= ' -j ' . LoadConfig::getParam($rH_cfg, 'pbdagcon', 'num_threads');
    $cmd .= ' ' . $infile;
    $cmd .= ' > ' . $outfile;
    $cmd .= ' && ';
    $cmd .= ' cat ' . $outfile . ' |'; # Then, convert to fastq.
    $cmd .= ' awk \'{if (\$0~/>/) {sub(/>/,\"@\",\$0);print;} else {l=length(\$0);q=\"\"; while (l--) {q=q \"9\"}printf(\"%s\n+\n%s\n\",\$0,q)}}\'';
    $cmd .= ' > ' . $outfileFastq;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub quiver {
  my $rH_cfg          = shift;
  my $infile          = shift;
  my $ref             = shift;
  my $outfileFasta    = shift;
  my $outfileVariants = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfileFasta]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh && ';
    $cmd .= ' memtime';
    $cmd .= ' quiver';
    $cmd .= ' -j ' . LoadConfig::getParam($rH_cfg, 'quiver', 'num_threads');
    $cmd .= ' ' . $infile;
    $cmd .= ' -r ' . $ref;
    $cmd .= ' -o ' . $outfileVariants;
    $cmd .= ' -o ' . $outfileFasta;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub callVariants {
  my $rH_cfg          = shift;
  my $infile          = shift;
  my $ref             = shift;
  my $outfileFasta    = shift;
  my $outfileVariants = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfileFasta, $outfileVariants]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' variantCaller.py';
    $cmd .= ' -j ' . LoadConfig::getParam($rH_cfg, 'quiver', 'num_threads');
    $cmd .= ' --algorithm quiver';
    $cmd .= ' --referenceFilename ' . $ref;
    $cmd .= ' --parameters best';
    $cmd .= ' -o ' . $outfileVariants;
    $cmd .= ' -o ' . $outfileFasta;
    $cmd .= ' ' . $infile;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub samtoh5 {
  my $rH_cfg          = shift;
  my $infile          = shift;
  my $ref             = shift;
  my $outfile         = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$infile], [$outfile]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' samtoh5';
    $cmd .= ' ' . $infile;
    $cmd .= ' ' . $ref;
    $cmd .= ' ' . $outfile;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub compareSequences {
  my $rH_cfg             = shift;
  my $cmpH5              = shift;
  my $controlRegionsFofn = shift;
  my $inputFofn          = shift;
  my $refUpload          = shift;
  my $tmpdir             = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$inputFofn], [$cmpH5]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' compareSequences.py';
    $cmd .= ' --info --useGuidedAlign --algorithm=blasr';
    $cmd .= ' --nproc=' . LoadConfig::getParam($rH_cfg, 'compareSequences', 'num_threads');
    $cmd .= ' --noXML --h5mode=w';
    $cmd .= ' --h5fn=' . $cmpH5;
    $cmd .= ' --seed=1 --minAccuracy=0.75 --minLength=50 --useQuality -x -minMatch 12 -x -bestn 10 -x -minPctIdentity 70.0 --placeRepeatsRandomly';
    $cmd .= ' --tmpDir=' . $tmpdir;
    $cmd .= ' --debug';
    $cmd .= ' --regionTable=' . $controlRegionsFofn;
    $cmd .= ' \"' . $inputFofn . '\"';
    $cmd .= ' \"' . $refUpload . '\"';
    $cmd .= ' > /dev/null';
    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub loadPulses {
  my $rH_cfg             = shift;
  my $cmpH5              = shift;
  my $inputFofn          = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$cmpH5, $inputFofn], [$cmpH5 . ".loadedPulses"]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' loadPulses';
    $cmd .= ' ' . $inputFofn;
    $cmd .= ' ' . $cmpH5;
    $cmd .= ' -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread';
    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub variantCaller {
  my $rH_cfg             = shift;
  my $cmpH5              = shift;
  my $refFasta           = shift;
  my $outfileVariants    = shift;
  my $outfileFasta       = shift;
  my $outfileFastq       = shift;

  my $outfileFastaUncompressed = $outfileFasta;
  $outfileFastaUncompressed =~ s/\.gz|\.gzip//;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$cmpH5 . ".loadedPulses"], [$outfileVariants, $outfileFasta, $outfileFastq]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' variantCaller.py';
    $cmd .= ' -P' . LoadConfig::getParam($rH_cfg, 'variantCaller', 'protocol');
    $cmd .= ' -v';
    $cmd .= ' -j' . LoadConfig::getParam($rH_cfg, 'variantCaller', 'num_threads');
    $cmd .= ' --algorithm=' . LoadConfig::getParam($rH_cfg, 'variantCaller', 'algorithm');;
    $cmd .= ' ' . $cmpH5;
    $cmd .= ' -r ' . $refFasta;
    $cmd .= ' -o ' . $outfileVariants;
    $cmd .= ' -o ' . $outfileFasta;
    $cmd .= ' -o ' . $outfileFastq;
    $cmd .= ' > /dev/null';
    $cmd .= ' && ';
    $cmd .= ' gunzip -c ' . $outfileFasta . ' > ' . $outfileFastaUncompressed;
    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub sortH5 {
  my $rH_cfg             = shift;
  my $cmpH5              = shift;
  my $cmpH5Out           = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$cmpH5], [$cmpH5Out]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' cmph5tools.py';
    $cmd .= ' -vv sort --deep --inPlace --outFile ' . $cmpH5Out;
    $cmd .= ' ' . $cmpH5;
    $cmd .= ' > /dev/null';
    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub runCASpecWriter {
  my $rH_cfg             = shift;
  my $genomeSize         = shift;
  my $estimatedCoverage  = shift;
  my $merSize            = shift;
  my $readsFasta         = shift;
  my $specOut            = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs(undef, undef);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' &&';
    $cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' runCASpecWriter.py';
    $cmd .= ' -vv ';
    $cmd .= ' --bitTable=\${SEYMOUR_HOME}/analysis/etc/celeraAssembler/bitTable';
    $cmd .= ' --interactiveTmpl=\${SEYMOUR_HOME}/analysis/etc/cluster/SGE/interactive.tmpl';
    $cmd .= ' --smrtpipeRc=\${SEYMOUR_HOME}/analysis/etc/smrtpipe.rc';
    $cmd .= ' --genomeSize ' . $genomeSize;
    $cmd .= ' --defaultFrgMinLen=' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'minReadSize');
    $cmd .= ' --xCoverage=' . $estimatedCoverage;
    $cmd .= ' --ovlErrorRate=' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'ovlErrorRate');
    $cmd .= ' --ovlMinLen=' . LoadConfig::getParam($rH_cfg, 'celeraConfig', 'ovlMinLen');
    $cmd .= ' --corrReadsFasta ' . $readsFasta;
    $cmd .= ' --specOut=' . $specOut;
    $cmd .= ' --merSize=' . $merSize;
    $cmd .= ' --sgeName=pacbioReads';
    $cmd .= ' --gridParams=\"useGrid:0,scriptOnGrid:0,frgCorrOnGrid:0,ovlCorrOnGrid:0\"';
    $cmd .= ' --maxSlotPerc=1';
    $cmd .= ' ' . LoadConfig::getParam($rH_cfg, 'default', 'celeraSettings');
    #$cmd .= ' \${SEYMOUR_HOME}/analysis/etc/celeraAssembler/template.spec';

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}

sub summarizePolishing {
  my $rH_cfg             = shift;
  my $sampleName         = shift;
  my $reference          = shift;
  my $alignedReadsCmpH5  = shift;
  my $alignmentSummary   = shift;
  my $coverageBed        = shift;
  my $alignedReadsSam    = shift;
  my $variantsGff        = shift;
  my $variantsBed        = shift;
  my $variantsVcf        = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$alignedReadsCmpH5, $reference], [$variantsVcf]);

  if (!$ro_job->isUp2Date()) {
    my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['smrtanalysis', 'moduleVersion.smrtanalysis']
    ]) . ' && source \${SEYMOUR_HOME}/etc/setup.sh &&';
    $cmd .= ' memtime';
    $cmd .= ' summarizeCoverage.py --reference ' . $reference . ' --numRegions=500 ' . $alignedReadsCmpH5 . ' > ' . $alignmentSummary;
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' gffToBed.py --name=meanCoverage --description=\"Mean coverage of genome in fixed interval regions\" coverage ' . $alignmentSummary . ' > ' . $coverageBed;
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' loadSequencingChemistryIntoCmpH5.py --xml ' . LoadConfig::getParam($rH_cfg, 'summarizePolishing', 'chemistryMapping') . ' --h5 ' . $alignedReadsCmpH5;
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' h5repack -f GZIP=1 ' . $alignedReadsCmpH5 . ' ' . $alignedReadsCmpH5 . '.repacked && mv ' . $alignedReadsCmpH5 . '.repacked ' . $alignedReadsCmpH5;
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' pbsamtools.py --bam --outfile ' . $alignedReadsSam . ' --refrepos ' . $reference . ' --readGroup movie ' . $alignedReadsCmpH5;
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' cmph5tools.py -vv sort --deep --inPlace ' . $alignedReadsCmpH5;
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' summarizeConsensus.py --variantsGff ' . $variantsGff . ' ' .  $alignmentSummary . '  --output ' . $alignmentSummary  . '.tmp && mv ' . $alignmentSummary . '.tmp ' . $alignmentSummary;
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' gffToBed.py --name=variants --description=\'PacBio: snps, insertions, and deletions derived from consensus calls against reference\' variants ' . $variantsGff . ' > ' . $variantsBed;
    $cmd .= ' && ';
    $cmd .= 'memtime';
    $cmd .= ' gffToVcf.py --globalReference=' . $sampleName . ' ' . $variantsGff . ' > ' . $variantsVcf;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
}
#  summarizeCoverage.py --reference /lb/scratch/jtrembla/pacbio_test/H41E5/101X/assembly/H41E5101X --numRegions=500 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/alignment_summary.gff || exit $?; echo "Task 0 completed at `date -u`" || exit $?;

#  gffToBed.py --name=meanCoverage --description="Mean coverage of genome in fixed interval regions" coverage /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/alignment_summary.gff > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/coverage.bed || exit $?; echo "Task 0 completed at `date -u`" || exit $?;

#  loadSequencingChemistryIntoCmpH5.py --xml /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/chemistry_mapping.xml --h5 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 || exit $?; echo "Task 0 completed at `date -u`" || exit $?;echo 'Chemistry Load Complete' || exit $?; echo "Task 1 completed at `date -u`" || exit $?;

#  ((which h5repack && (h5repack -f GZIP=1 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5_TMP && mv /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5_TMP /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5)) || echo 'no h5repack found, continuing w/out') || exit $?; echo "Task 0 completed at `date -u`" || exit $?;

#  pbsamtools.py --bam --outfile /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.sam --refrepos /lb/scratch/jtrembla/pacbio_test/H41E5/101X/assembly/H41E5101X --readGroup movie /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 || exit $?; echo "Task 0 completed at `date -u`" || exit $?;

#  cmph5tools.py -vv sort --deep --inPlace /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 || exit $?; echo "Task 0 completed at `date -u`" || exit $?;echo 'Sorting Complete' || exit $?; echo "Task 1 completed at `date -u`" || exit $?;

#  extractUnmappedSubreads.py /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/filtered_subreads.fasta /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/aligned_reads.cmp.h5 /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/control_reads.cmp.h5 > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/unmappedSubreads.fasta || exit $?; echo "Task 0 completed at `date -u`" || exit $?;

#  summarizeConsensus.py --variantsGff /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/alignment_summary.gff --output /lb/scratch/jtrembla/tmp/tmppIn9qP.gff || exit $?; echo "Task 0 completed at `date -u`" || exit $?;mv /lb/scratch/jtrembla/tmp/tmppIn9qP.gff /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/alignment_summary.gff || exit $?; echo "Task 1 completed at `date -u`" || exit $?;

#  gffToBed.py --name=variants --description='PacBio: snps, insertions, and deletions derived from consensus calls against reference' variants /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.bed || exit $?; echo "Task 0 completed at `date -u`" || exit $?;

#  gffToVcf.py --globalReference=H41E5101X /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff > /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.vcf || exit $?; echo "Task 0 completed at `date -u`" || exit $?;

#  gzip -f /lb/scratch/jtrembla/pacbio_test/H41E5/101X/polishing/data/variants.gff || exit $?; echo "Task 0 completed at `date -u`" || exit $?;
  #########

1;
