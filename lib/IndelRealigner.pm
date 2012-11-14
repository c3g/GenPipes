#!/usr/env/perl

=head1 NAME

I<Trimmomatic>

=head1 SYNOPSIS

Trimmomatic->trim()

=head1 DESCRIPTION

B<Trimmomatic> is a library that trims fastqs

Input = file_name

Output = array


=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package IndelRealigner;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------
sub realign {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $seqName         = shift;
  my $processUnmapped = shift;

  print "mkdir -p $sampleName/realign\n";
  print "JOB_IDS=\"\"\n";

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'referenceFasta');
  my $sortedBAM = $sampleName.'/'.$sampleName.'.sorted.bam';
  my $intervalOutput = $sampleName.'/realign/'.$seqName.'.intervals';
  my $realignOutput = $sampleName.'/realign/'.$seqName.'.bam';
  
  my $intervalCommand;
  $intervalCommand .= 'module load mugqic/GenomeAnalysisTKLite/2.1-13 ;';
  $intervalCommand .= ' java -Xmx1G -jar ${GATK_JAR}';
  $intervalCommand .= ' -T RealignerTargetCreator';
  $intervalCommand .= ' -R '.$refGenome;
  $intervalCommand .= ' -o '.$intervalOutput;
  $intervalCommand .= ' -I '.$sortedBAM;
  $intervalCommand .= ' -L '.$seqName;

  my $realignCommand;
  $realignCommand .= 'module load mugqic/GenomeAnalysisTKLite/2.1-13 ;';
  $realignCommand .= ' java -Xmx1G -jar ${GATK_JAR}';
  $realignCommand .= ' -T IndelRealigner';
  $realignCommand .= ' -R '.$refGenome;
  $realignCommand .= ' -targetIntervals '.$intervalOutput;
  $realignCommand .= ' -o '.$intervalOutput;
  $realignCommand .= ' -I '.$sortedBAM;
  $realignCommand .= ' -L '.$seqName;
  if($processUnmapped == 1) {
    $realignCommand .= ' -L unmapped';
  }

    my $TARGET_MOAB_COMMAND="msub -d $CURRENT_DIR -V -l walltime=24:00:0 -l nodes=1:ppn=2 -j oe -q sw -o $output_jobs $dependency -m ae -M $EMAIL_NOTIFICATION -N int$sampleName"."_"."$chr";
    if($c == 0) {
       -L  -o $sampleName/realign/$chr.bam";
      $REALIGNER_MOAB_COMMAND="msub -d $CURRENT_DIR -V -l walltime=48:00:0 -l nodes=1:ppn=6 -j oe -o $output_jobs -m ae -M $EMAIL_NOTIFICATION -q sw -N realn$sampleName"."_"."$chr -W x=depend:afterok:\${INTERVAL_JOB_ID}";
    }
    else {
      $REALIGNER_MOAB_COMMAND="msub -d $CURRENT_DIR -V -l walltime=48:00:0 -l nodes=1:ppn=6 -j oe -o $output_jobs -m ae -M $EMAIL_NOTIFICATION -q sw -N realn$sampleName"."_"."$chr -W x=depend:afterok:\${INTERVAL_JOB_ID}";
    }
    #my $REALIGNER_MOAB_COMMAND="msub -d $CURRENT_DIR -V -l walltime=48:00:0 -l nodes=1:ppn=6 -j oe -o $output_jobs -m ae -M $EMAIL_NOTIFICATION -q sw -N realn$sampleName"."_"."$chr -W x=depend:afterok:\${INTERVAL_JOB_ID}";
    unless (-e "$sampleName/realign/$chr.intervals.done" && -e "$sampleName/realign/$chr.bam.done") {
      print "INTERVAL_JOB_ID=`echo \"rm -f $sampleName/realign/$chr.intervals.done ; $JAVA_BIN -Djava.io.tmpdir=$TMP_DIR -Xmx4G -jar $GATK_JAR $TARGET_OPTIONS && touch $sampleName/realign/$chr.intervals.done\" | $TARGET_MOAB_COMMAND | grep \"[0-9]\"`\n";
      print "JOB_IDS=\${JOB_IDS}:`echo \"rm -f $sampleName/realign/$chr.bam.done ; $JAVA_BIN -Djava.io.tmpdir=$TMP_DIR -Xmx11G -jar $GATK_JAR $REALIGNER_OPTIONS && touch $sampleName/realign/$chr.bam.done\" | $REALIGNER_MOAB_COMMAND | grep \"[0-9]\"`\n";
      print "FILES_TO_MERGE=\"\${FILES_TO_MERGE} INPUT=$sampleName/realign/$chr.bam\"\n";
    }
  }














  my $latestBam;
  my $bamInputs;
  my $countInputs;
  my $outputBAM = $sampleName.'/'.$sampleName.'.sorted.bam';
  for my $rH_laneInfo (@$rAoH_sampleLanes) {
    my $directory = $sampleName."/run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'}."/";
    my $sortedLaneBamFile = $directory.$rH_laneInfo->{'name'}.".sorted.bam";
    my $laneStatsFile = $directory.$rH_laneInfo->{'name'}.".counts";
    my $runName = $sampleName."_run".$rH_laneInfo->{'runId'}."_".$rH_laneInfo->{'lane'};


    if(!defined($latestBam)) {
      $latestBam = -M $sortedLaneBamFile;
    }
    else {
      my $modDate = -M $sortedLaneBamFile;
      if($modDate > $latestBam) {
        $latestBam = $modDate;
      }
    }
    $bamInputs .= 'INPUT='.$sortedLaneBamFile.' ';
    $countInputs .= $laneStatsFile.' ';
  }

  my $command;
  if(!defined($latestBam) || !defined(-M $outputBAM) || $latestBam > -M $outputBAM) {
    $command .= 'module load mugqic/picard/1.77 ;';
    $command .= ' java -Xmx'.LoadConfig::getParam($rH_cfg, 'mergeLanes', 'mergeRam').' -jar ${PICARD_HOME}/MergeSamFiles.jar';
    $command .= ' VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true';
    $command .= ' '.$bamInputs;
    $command .= ' OUTPUT='.$outputBAM;
    $command .= ' MAX_RECORDS_IN_RAM='.LoadConfig::getParam($rH_cfg, 'mergeLanes', 'mergeRecInRam');
    $command .= ' ; cat '.$countInputs.' > $sampleName/$sampleName.runLane.counts';
  }
  return $command;
}

1;
