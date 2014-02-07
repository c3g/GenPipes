#!/usr/env/perl

=head1 NAME

I<BWA>

=head1 SYNOPSIS

BWA->aln()

=head1 DESCRIPTION

B<BWA> is a library that aligns fastqs on a reference genome

Input = file_name

Output = array


=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package BWA;

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
sub mem {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $pair1           = shift;
  my $pair2           = shift;
  my $single          = shift;
  my $optOutputPrefix = shift;
  my $rgId            = shift;
  my $rgSample        = shift;
  my $rgLibrary       = shift;
  my $rgPlatformUnit  = shift;
  my $rgCenter        = shift;
  my $bwaRefIndex     = shift;

  if (!defined($bwaRefIndex)) {
    $bwaRefIndex = LoadConfig::getParam( $rH_cfg, 'mem', 'bwaRefIndex' );
  }
    
  my $outputBAM = $optOutputPrefix.'.sorted.bam';

  my $rA_inputs;
  my $dateToTest;
  if(defined($pair1) && defined($pair2)) {
    $rA_inputs = [$pair1, $pair2];
  }
  else {
    $rA_inputs = [$single];
  }

  my $ro_job = new Job();
  $ro_job->testInputOutputs($rA_inputs, [$outputBAM]);

  if (!$ro_job->isUp2Date()) {
    my $rgTag = "'" . '@RG\tID:' . $rgId . '\tSM:' . $rgSample . '\tLB:' . $rgLibrary . '\tPU:run' . $rgPlatformUnit . '\tCN:' . $rgCenter . '\tPL:Illumina' . "'";
    my $bwaCommand;
    $bwaCommand .= 'module load '.LoadConfig::getParam($rH_cfg, 'mem', 'moduleVersion.bwa').' '.LoadConfig::getParam($rH_cfg, 'mem', 'moduleVersion.picard').' '.LoadConfig::getParam($rH_cfg, 'mem', 'moduleVersion.java').' &&';
    $bwaCommand .= ' bwa mem ';
    $bwaCommand .= ' '.LoadConfig::getParam( $rH_cfg, 'mem', 'bwaExtraFlags' );
    $bwaCommand .= ' -R ' . $rgTag;
    $bwaCommand .= ' ' . $bwaRefIndex;
    if ( defined($pair1) && defined($pair2) ) {
      $bwaCommand .= ' ' . $pair1;
      $bwaCommand .= ' ' . $pair2;
    }
    elsif ( defined($single) ) {
      $bwaCommand .= ' ' . $single;
    }
    else {
      die "Unknown runType, not paired or single\n";
    }
    $bwaCommand .= ' | java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'mem', 'tmpDir');
    $bwaCommand .= ' '.LoadConfig::getParam($rH_cfg, 'mem', 'extraJavaFlags');
    $bwaCommand .= ' -Xmx'.LoadConfig::getParam($rH_cfg, 'mem', 'sortRam');
    $bwaCommand .= ' -jar \${PICARD_HOME}/SortSam.jar';
    $bwaCommand .= '  INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate';
    $bwaCommand .= ' OUTPUT=' . $outputBAM;
    $bwaCommand .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam( $rH_cfg, 'mem', 'sortRecInRam' );

    $ro_job->addCommand($bwaCommand);
  }

  return $ro_job;
}

sub aln {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $pair1           = shift;
  my $pair2           = shift;
  my $single          = shift;
  my $optOutputPrefix = shift;
  my $rgId            = shift;
  my $rgSample        = shift;
  my $rgLibrary       = shift;
  my $rgPlatformUnit  = shift;
  my $rgCenter        = shift;
  my $indexToUse      = shift;

  my $rO_job;
  if ( defined($pair1) && defined($pair2) ) {
    $rO_job = pairCommand( $rH_cfg, $sampleName, $pair1, $pair2, $optOutputPrefix, $rgId, $rgSample, $rgLibrary, $rgPlatformUnit, $rgCenter, $indexToUse);
  }
  elsif ( defined($single) ) {
    $rO_job = singleCommand( $rH_cfg, $sampleName, $single, $optOutputPrefix, $rgId, $rgSample, $rgLibrary, $rgPlatformUnit, $rgCenter, $indexToUse);
  }
  else {
    die "Unknown runType, not paired or single\n";
  }

  return $rO_job;
}

sub pairCommand {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $pair1           = shift;
  my $pair2           = shift;
  my $optOutputPrefix = shift;
  my $rgId            = shift;
  my $rgSample        = shift;
  my $rgLibrary       = shift;
  my $rgPlatformUnit  = shift;
  my $rgCenter        = shift;
  my $indexToUse      = shift;

  my $bwaRefIndex = LoadConfig::getParam( $rH_cfg, 'aln', 'bwaRefIndex' );
  if (defined $indexToUse) {
    $bwaRefIndex = $indexToUse;
  }
  
  my $outputSai1Name = $optOutputPrefix.'.pair1.sai';
  my $outputSai2Name = $optOutputPrefix.'.pair2.sai';
  my $outputBAM = $optOutputPrefix.'.sorted.bam';

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$pair1, $pair2], [$outputSai1Name, $outputSai2Name, $outputBAM]);

  if (!$ro_job->isUp2Date()) {
    my $sai1Command = "";
    my $sai2Command = "";
    my $bwaCommand  = "";
    $sai1Command .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' &&';
    $sai1Command .= 'bwa aln';
    $sai1Command .= ' -t ' . LoadConfig::getParam( $rH_cfg, 'aln', 'bwaAlnThreads' );
    $sai1Command .= ' ' . $bwaRefIndex;
    $sai1Command .= ' ' . $pair1;
    $sai1Command .= ' -f ' . $outputSai1Name;
    $ro_job->addCommand($sai1Command);

    $sai2Command .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' &&';
    $sai2Command .= 'bwa aln';
    $sai2Command .= ' -t ' . LoadConfig::getParam( $rH_cfg, 'aln', 'bwaAlnThreads' );
    $sai2Command .= ' ' . $bwaRefIndex;
    $sai2Command .= ' ' . $pair2;
    $sai2Command .= ' -f ' . $outputSai2Name;
    $ro_job->addCommand($sai2Command);

    my $rgTag = "'" . '@RG\tID:' . $rgId . '\tSM:' . $rgSample . '\tLB:' . $rgLibrary . '\tPU:run' . $rgPlatformUnit . '\tCN:' . $rgCenter . '\tPL:Illumina' . "'";
    $bwaCommand .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.picard').' '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.java').' &&';
    $bwaCommand .= ' bwa sampe ';
    $bwaCommand .= ' '.LoadConfig::getParam( $rH_cfg, 'aln', 'bwaExtraSamXeFlags' );
    $bwaCommand .= ' -r ' . $rgTag;
    $bwaCommand .= ' ' . $bwaRefIndex;
    $bwaCommand .= ' ' . $outputSai1Name;
    $bwaCommand .= ' ' . $outputSai2Name;
    $bwaCommand .= ' ' . $pair1;
    $bwaCommand .= ' ' . $pair2;
    $bwaCommand .= ' | java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'aln', 'tmpDir');
    $bwaCommand .= ' '.LoadConfig::getParam($rH_cfg, 'aln', 'extraJavaFlags');
    $bwaCommand .= ' -Xmx'.LoadConfig::getParam($rH_cfg, 'aln', 'sortRam');
    $bwaCommand .= ' -jar \${PICARD_HOME}/SortSam.jar';
    $bwaCommand .= '  INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate';
    $bwaCommand .= ' OUTPUT=' . $outputBAM;
    $bwaCommand .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam( $rH_cfg, 'aln', 'sortRecInRam' );

    $ro_job->addCommand($bwaCommand);
  }

  return $ro_job;
}

sub singleCommand {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $single          = shift;
  my $optOutputPrefix = shift;
  my $rgId            = shift;
  my $rgSample        = shift;
  my $rgLibrary       = shift;
  my $rgPlatformUnit  = shift;
  my $rgCenter        = shift;
  my $indexToUse      = shift;
    
  my $bwaRefIndex = LoadConfig::getParam( $rH_cfg, 'aln', 'bwaRefIndex' );
  if (defined $indexToUse){
    $bwaRefIndex = $indexToUse;
  }
  
  my $outputSaiName = $optOutputPrefix.'.single.sai';
  my $outputBAM = $optOutputPrefix.'.sorted.bam';

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$single], [$outputSaiName, $outputBAM]);

  # we could test sai and bam separately...
  if (!$ro_job->isUp2Date()) {
    my $saiCommand = "";

    $saiCommand .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' &&';
    $saiCommand .= 'bwa aln';
    $saiCommand .= ' -t ' . LoadConfig::getParam( $rH_cfg, 'aln', 'bwaAlnThreads' );
    $saiCommand .= ' ' . $bwaRefIndex;
    $saiCommand .= ' ' . $single;
    $saiCommand .= ' -f ' . $outputSaiName;
    $ro_job->addCommand($saiCommand);

    my $rgTag = "'" . '@RG\tID:' . $rgId . '\tSM:' . $rgSample . '\tLB:' . $rgLibrary . '\tPU:run' . $rgPlatformUnit . '\tCN:' . $rgCenter . '\tPL:Illumina' . "'";
    my $bwaCommand = "";
    $bwaCommand .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.picard').' '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.java').' &&';
    $bwaCommand .= ' bwa samse';
    $bwaCommand .= ' '.LoadConfig::getParam( $rH_cfg, 'aln', 'bwaExtraSamXeFlags' );
    $bwaCommand .= ' -r ' . $rgTag;
    $bwaCommand .= ' ' . $bwaRefIndex;
    $bwaCommand .= ' ' . $outputSaiName;
    $bwaCommand .= ' ' . $single;
    $bwaCommand .= ' | java -Djava.io.tmpdir='.LoadConfig::getParam($rH_cfg, 'aln', 'tmpDir');
    $bwaCommand .= ' '.LoadConfig::getParam($rH_cfg, 'aln', 'extraJavaFlags');
    $bwaCommand .= ' -Xmx'.LoadConfig::getParam($rH_cfg, 'aln', 'sortRam');
    $bwaCommand .= ' -jar \${PICARD_HOME}/SortSam.jar';
    $bwaCommand .= ' INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate';
    $bwaCommand .= ' OUTPUT=' . $outputBAM;
    $bwaCommand .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam( $rH_cfg, 'aln', 'sortRecInRam' );

    $ro_job->addCommand($bwaCommand);
  }

  return $ro_job;
}

sub index {
  my $rH_cfg   = shift;
  my $toIndex  = shift;

  my $ro_job = new Job();
  $ro_job->testInputOutputs([$toIndex], [$toIndex.'.bwt']);

  if (!$ro_job->isUp2Date()) {
    my $command = 'module load '.LoadConfig::getParam($rH_cfg, 'index', 'moduleVersion.bwa').' &&';
    $command .= ' bwa index ' . $toIndex;

    $ro_job->addCommand($command);
  }
  return $ro_job
}

1;

