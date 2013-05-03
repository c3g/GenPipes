#!/usr/bin/perl

require 5.008;

use strict;
use POSIX qw( strftime );
use File::Copy;
use File::Find;
use XML::Simple;
use Data::Dumper;

use Getopt::Long;

my $version = "0.2";

&main();

sub main {
  my $runDirectory;
  my @to;
  my $verbose;
  my $email;
  my $userMask;
  my $baseCallDir;
  my $fromEmail;
  my $nanuqAuthFile;
  my $isMiSeq = 0;

  my $result = GetOptions(
    "from=s"           => \$fromEmail,
    "to=s"             => \@to,
    "runDir=s"         => \$runDirectory,
    "nanuqAuthFilei=s" => \$nanuqAuthFile,
    "mask=s"           => \$userMask,
    "email=s"          => \$email,
    "baseCallDir=s"    => \$baseCallDir,
    "verbose"          => \$verbose
  );

  my $errMsg = "";
  if(!defined($nanuqAuthFile) || !-e $nanuqAuthFile) {
    $errMsg .= "Missing nanuqAuthFile\n";
  }
  if(!defined($email) || length($email) == 0) {
    $errMsg .= "Missing email\n";
  }
  if(length($errMsg)) {
    die $errMsg;
  }


  if(!defined($baseCallDir)) {
    $baseCallDir = $runDirectory.'/Data/Intensities/BaseCalls';
  }

  my $runID;
  my $runName;
  if ( $runDirectory =~ /.*_\d+HS\d\d[AB]/ ) {
    ($runName,$runID) = $runDirectory =~ /.*\/(\d+_[^_]+_\d+_[^_]+_(\d+)HS.+)/;
  }
  elsif ( $runDirectory =~ /.*\d+_[^_]+_\d+_.+/ ) {
    ($runName,$runID) = $runDirectory =~ /.*\/(\d+_([^_]+_\d+)_.*)/;
  }
  else {
    die "Not a valid run folder: " . $runDirectory . "\n";
  }

  if ( !-e $runDirectory . "/RunInfo.xml" ) {
    die "Missing RunInfo.xml: " . $runDirectory . "\n";
  }

  my $xml     = new XML::Simple( 'ForceArray' => ['Read'] );
  my $runInfo = $xml->XMLin( $runDirectory . "/RunInfo.xml" );
  my $nbLanes = $runInfo->{'Run'}->{'FlowcellLayout'}->{'LaneCount'};
  my $nbReads = @{ $runInfo->{'Run'}->{'Reads'}->{'Read'} };
  my $fileToLookFor;
  my $nbTemplateReads = 0;
  my @parsedReads;
  my $mask;
  for my $rh_read (@{$runInfo->{'Run'}->{'Reads'}->{'Read'}}) {
    my $nbCycles = $rh_read->{'NumCycles'};
    my $isIndexed = $rh_read->{'IsIndexedRead'};
    push(@parsedReads, {nbCycles=>$nbCycles , isIndexed=>$isIndexed});
    if(length($mask) != 0) {
      $mask .= ',';
    }
    if($isIndexed eq "Y") {
      $mask .= 'I'.$nbCycles;
    }
    else {
      $mask .= 'Y'.$nbCycles;
      $nbTemplateReads++;
    }
  }

  my @readsInfo;
  if ( defined($userMask) ) {
    my @values = split(',', $userMask);
    for my $value (@values) {
      if($value =~ /^[yY](\d+)/) {
        push(@readsInfo, {nbCycles=>$1 , isIndexed=>'N'});
      }
      elsif ($value =~ /^[I](\d+)/) {
        push(@readsInfo, {nbCycles=>$1 , isIndexed=>'Y'});
      }
      else {
        die "Unknown mask character: ".$value."\n";
      }
    }
    $mask = $userMask;
  }
  else {
    @readsInfo = @parsedReads;
  }

  my $indexLength = computeIndexLength(\@readsInfo);
  
  if($runDirectory =~ /_M00/){
    $isMiSeq = 1;
  }
  my $rHoH_samples = getSampleSheet( $email, $isMiSeq, $runDirectory, $nanuqAuthFile, $runID, \@readsInfo, $indexLength );
  generateIndexCounts($email, $runDirectory, $baseCallDir, $runID, $runName, $nbLanes, $indexLength, \@readsInfo);
  generateFastq($email, $runDirectory, $baseCallDir, $runID, $mask);
  generateMD5($email, $runDirectory, $baseCallDir, $runID, $nbTemplateReads, $rHoH_samples);
  generateQCGraphs($email, $runDirectory, $baseCallDir, $runID, $nbTemplateReads, $rHoH_samples);
  generateBlasts($email, $runDirectory, $baseCallDir, $runID, $nbTemplateReads, $rHoH_samples);
  
  my $destinationFolder = '/lb/robot/hiSeqSequencer/hiSeqRuns-drop/';
  if($isMiSeq == 1) {
    $destinationFolder = '/lb/robot/miSeqSequencer/miSeqRuns-drop/';
  }

  
  print 'echo "rsync -avP --include \"**/*onfig*\" --include \"Unaligned/**\" --exclude \"Thumbnail_Images/\" --exclude \"Data/Intensities/B*/*\" --include \"Data/Intensities/B*/\" --exclude \"Data/Intensities/*\" '.$runDirectory.' '.$destinationFolder.' ; setfacl -R -m mask:rwx '.$destinationFolder.$runName.'" | msub -l qos=hiseq -d '.$runDirectory.' -V -l nodes=1:ppn=1 -l walltime=24:00:0 -q sw -j oe -N rsync.'.$runID.' -W x=depend:afterok:${QC_JOB_IDS}:${BARCODE_ID_JOB_IDS}:${BLAST_JOB_IDS} -m ae -M '.$email."\n";
#  print 'echo "echo \"RUN='.$runID.' ; rsync --delete -avP -e \"ssh -c arcfour\" --include \"**/*onfig*\" --exclude \"Unaligned/**.old\" --include \"Unaligned/**\" --exclude \"Thumbnail_Images/\" --exclude \"Data/Intensities/B*/*\" --include \"Data/Intensities/B*/\" --exclude \"Data/Intensities/*\" abacus:\${RUN} /data/newrobot/hiSeqSequencer/hiSeqRuns-drop/ ; chmod -R ug+rw \${RUN} ; find \${RUN} -type d | xargs chmod ug+x ; time setfacl -R -m mask:rwx \${RUN} ; echo \"\${RUN}\"" | msub -l qos=hiseq -d '.$runDirectory.' -V -l nodes=1:ppn=1 -l walltime=1:00:0 -q sw -j oe -N final.'.$runID.' -W x=depend:afterok:${QC_JOB_IDS}:${BARCODE_ID_JOB_IDS}:${BLAST_JOB_IDS} -m ae -M '.$email."\n";
}

sub generateBlasts {
  my $email          = shift;
  my $runDirectory = shift;
  my $baseCallDir = shift;
  my $runID = shift;
  my $nbTemplateReads = shift;
  my $rHoH_samples = shift;
  
  print "BLAST_JOB_IDS=\"\"\n";
  for my $rH_sample (values(%{$rHoH_samples})) {
    if($nbTemplateReads == 1) {
      print 'BLAST_JOB_IDS=${BLAST_JOB_IDS}:`echo "mkdir -p '.$runDirectory.'/Unaligned/Blast_sample ; runBlast.sh 100000 '.$runDirectory.'/Unaligned/Blast_sample/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.' '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R1_001.fastq.gz" | msub -l qos=hiseq -d '.$runDirectory.'/Unaligned -V -l nodes=1:ppn=2 -l walltime=6:00:0 -q sw -j oe -N blast.'.$runID.'.'.$rH_sample->{'name'}.' -W x=depend:afterok:${FASTQ_JOB_ID} -m ae -M '.$email.' | grep "[0-9]"`'."\n";
    } else {
      print 'BLAST_JOB_IDS=${BLAST_JOB_IDS}:`echo "mkdir -p '.$runDirectory.'/Unaligned/Blast_sample ; runBlast.sh 100000 '.$runDirectory.'/Unaligned/Blast_sample/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.' '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R1_001.fastq.gz '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R2_001.fastq.gz" | msub -l qos=hiseq -d '.$runDirectory.'/Unaligned -V -l nodes=1:ppn=2 -l walltime=6:00:0 -q sw -j oe -N blast.'.$runID.'.'.$rH_sample->{'name'}.' -W x=depend:afterok:${FASTQ_JOB_ID} -m ae -M '.$email.' | grep "[0-9]"`'."\n";
    }
  }
}

sub generateMD5 {
  my $email          = shift;
  my $runDirectory = shift;
  my $baseCallDir = shift;
  my $runID = shift;
  my $nbTemplateReads = shift;
  my $rHoH_samples = shift;
  
  print "MD5_JOB_IDS=\"\"\n";
  for my $rH_sample (values(%{$rHoH_samples})) {
    if($nbTemplateReads == 1) {
      print 'MD5_JOB_IDS=${MD5_JOB_IDS}:`echo "md5sum -b '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R1_001.fastq.gz > '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R1_001.fastq.gz.md5" | msub -l qos=hiseq -d '.$runDirectory.'/Unaligned -V -l nodes=1:ppn=5 -l walltime=12:00:0 -q sw -j oe -N md5.'.$runID.'.'.$rH_sample->{'name'}.' -W x=depend:afterok:${FASTQ_JOB_ID} -m ae -M '.$email.' | grep "[0-9]"`'."\n";
    } else {
      print 'MD5_JOB_IDS=${MD5_JOB_IDS}:`echo "md5sum -b '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R1_001.fastq.gz > '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R1_001.fastq.gz.md5 ; md5sum -b '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R2_001.fastq.gz > '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R2_001.fastq.gz.md5" | msub -l qos=hiseq -d '.$runDirectory.'/Unaligned -V -l nodes=1:ppn=5 -l walltime=12:00:0 -q sw -j oe -N md5.'.$runID.'.'.$rH_sample->{'name'}.' -W x=depend:afterok:${FASTQ_JOB_ID} -m ae -M '.$email.' | grep "[0-9]"`'."\n";
    }
  }
}

sub generateQCGraphs {
  my $email          = shift;
  my $runDirectory = shift;
  my $baseCallDir = shift;
  my $runID = shift;
  my $nbTemplateReads = shift;
  my $rHoH_samples = shift;
  
  print "QC_JOB_IDS=\"\"\n";
  for my $rH_sample (values(%{$rHoH_samples})) {
    if($nbTemplateReads == 1) {
      print 'QC_JOB_IDS=${QC_JOB_IDS}:`echo "mkdir -p '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/qc ; time java -Xmx20G -Djava.awt.headless=true -jar /sb/programs/analyste/software/mps-tools.jar -i '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R1_001.fastq.gz -Q33 -t FASTQ -o '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/qc/ -n5 -outputIlmnHandler -regionName '.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'" | msub -l qos=hiseq -d '.$runDirectory.'/Unaligned -V -l nodes=1:ppn=5 -l walltime=12:00:0 -q sw -j oe -N qc.'.$runID.'.'.$rH_sample->{'name'}.' -W x=depend:afterok:${FASTQ_JOB_ID} -m ae -M '.$email.' | grep "[0-9]"`'."\n";
    } else {
      print 'QC_JOB_IDS=${QC_JOB_IDS}:`echo "mkdir -p '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/qc;time java -Xmx20G -Djava.awt.headless=true -jar /sb/programs/analyste/software/mps-tools.jar -i '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R1_001.fastq.gz -i2 '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/'.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'_R2_001.fastq.gz -Q33 -t FASTQ -o '.$runDirectory.'/Unaligned/Project_nanuq/Sample_'.$rH_sample->{'name'}.'/qc -n5 -outputIlmnHandler -regionName '.$rH_sample->{'name'}.'_'.$rH_sample->{'index'}.'_L00'.$rH_sample->{'lane'}.'" | msub -l qos=hiseq -d '.$runDirectory.'/Unaligned -V -l nodes=1:ppn=5 -l walltime=12:00:0 -q sw -j oe -N qc.'.$runID.'.'.$rH_sample->{'name'}.' -W x=depend:afterok:${FASTQ_JOB_ID} -m ae -M '.$email.' | grep "[0-9]"`'."\n";
    }
  }
}

sub findFiles {
  my $rA_fastqs = shift;
  my $file = $File::Find::name;

  return unless -f $file;             # process files (-f), not directories
  return unless $file =~ m/(.*Unaligned\/Project_nanuq\/.*\/.*_R1_001.fastq.gz)/o;
  push(@{$rA_fastqs}, $file);
}

sub generateFastq {
  my $email          = shift;
  my $runDirectory = shift;
  my $baseCallDir = shift;
  my $runID = shift;
  my $mask = shift;

  print "# Generate Unaligned directory\n";
  if($mask !~ /I/) {
    print "configureBclToFastq.pl --input-dir $baseCallDir --sample-sheet $runDirectory/nanuqSampleSheet.csv --fastq-cluster-count 0\n";
  } else {
    print "configureBclToFastq.pl --input-dir $baseCallDir --sample-sheet $runDirectory/nanuqSampleSheet.csv --fastq-cluster-count 0 --mismatches 1 --use-bases-mask $mask\n";
  }
  
  print 'FASTQ_JOB_ID=`echo "make -j 12" | msub -l qos=hiseq -d '.$runDirectory.'/Unaligned -V -l nodes=1:ppn=12 -l walltime=12:00:0 -q sw -j oe -N fastq'.$runID.' -m ae -M '.$email.' | grep "[0-9]"`'."\n";
}

sub computeIndexLength {
  my $rAoH_readsInfo = shift;
  
  my $indexLength = 0;
  for my $rH_readInfo (@{$rAoH_readsInfo}) {
    if($rH_readInfo->{'isIndexed'} eq 'Y') {
      $indexLength += $rH_readInfo->{'nbCycles'};
    }
  }
  return $indexLength;
}

sub generateIndexCounts {
  my $email           = shift;
  my $runDirectory    = shift;
  my $baseCallDir     = shift;
  my $runID           = shift;
  my $runName         = shift;
  my $nbLanes         = shift;
  my $indexLength     = shift;
  my $rAoH_readsInfo  = shift;

  my $mask = "";
  my $indexPrinted=0;
  for my $rH_readInfo (@{$rAoH_readsInfo}) {
    if($rH_readInfo->{'isIndexed'} eq 'Y') {
      if($indexPrinted == 0) {
        $mask .= $indexLength.'B';
        $indexPrinted=1;
      } 
    }
    else {
      $mask .= $rH_readInfo->{'nbCycles'}.'T';
    }
  }

  if( $indexPrinted == 0) {
    print "# No Indexes, *NOT* Generating index counts\n";
  } else {
    print "# Generating index counts\n";
    
    print 'BARCODE_ID_JOB_IDS=""'."\n";
    for(my $lane=1; $lane <= $nbLanes; $lane++) {
      print 'BARCODE_ID_JOB_IDS=${BARCODE_ID_JOB_IDS}:`echo "java -Xmx55G -jar /sb/programs/analyste/software/CountIlluminaBarcodes-0.2-alpha-jar-with-dependencies.jar BASECALLS_DIR='.$baseCallDir.' LANE='.$lane.' READ_STRUCTURE='.$mask.' METRICS_FILE='.$runDirectory.'/'.$runName.'_'.$lane.'.metrics MAX_MISMATCHES=1 NUM_PROCESSORS=12 BARCODE_FILE=/sb/programs/analyste/software/barcodes.txt" | msub -l qos=hiseq -d '.$runDirectory.' -V -l nodes=1:ppn=12 -l walltime=12:00:0 -q sw -j oe -N idx_'.$runID.'_'.$lane.' -m ae -M '.$email.' | grep "[0-9]"`'."\n";
    }
  }
}

sub getSampleSheet {
  my $email          = shift;
  my $isMiSeq        = shift;
  my $runDirectory   = shift;
  my $nanuqAuthFile  = shift;
  my $runID          = shift;
  my $rAoH_readsInfo = shift;
  my $indexLength    = shift;

  if ( !-e $runDirectory . "/SampleSheet.nanuq.csv" ) {
    my $instrument = 'HiSeq';
    if($isMiSeq == 1) {
      $instrument = 'MiSeq';
    }
    print '#wget --no-cookies --directory-prefix '.$runDirectory.'/ --post-file '.$nanuqAuthFile.' https://genomequebec.mcgill.ca/nanuqMPS/sampleSheet/'.$instrument.'/'.$runID.'/SampleSheet.nanuq.csv'."\n";
    system('wget --no-cookies --directory-prefix '.$runDirectory.'/ --post-file '.$nanuqAuthFile.' https://genomequebec.mcgill.ca/nanuqMPS/sampleSheet/'.$instrument.'/'.$runID.'/SampleSheet.nanuq.csv');
    if ($? == -1) {
      print "failed to execute: $!\n";
      exit(1);
    }
    elsif ($? & 127) {
      printf "child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without';
      exit(1);
    }
    else {
      my $childValue = $? >> 8;
      if($childValue != 0) {
        printf "child exited with value %d\n", $childValue;
        exit(1);
      }
    }
  }

  my %retVal;
  open( SSHEET, $runDirectory . "/SampleSheet.nanuq.csv" ) or die "Can't open sample sheet: " . $runDirectory . "/SampleSheet.nanuq.csv\n";
  open( NANUQ_SHEET, '>', $runDirectory."/nanuqSampleSheet.csv" ) or die "Can't write nanuq sample sheet: " . $runDirectory . "/nanuqSampleSheet.csv\n";
  my $line = <SSHEET>;
  print NANUQ_SHEET $line;
  my @values = split( /,/, $line );
  my $indexColumnIdx = -1;
  my $laneColumnIdx = -1;
  my $sampleIDColumnIdx = -1;

  for ( my $idx = 0 ; $idx < @values ; $idx++ ) {
    if ( $values[$idx] eq 'Index' ) {
      $indexColumnIdx = $idx;
    }
    elsif ( $values[$idx] eq 'Lane' ) {
      $laneColumnIdx = $idx;
    }
    elsif ( $values[$idx] eq 'SampleID' ) {
      $sampleIDColumnIdx = $idx;
    }
  }
  if ( $indexColumnIdx == -1 || $laneColumnIdx == -1 || $sampleIDColumnIdx == -1) {
    die "Missing columns\n";
  }

    while ( my $line = <SSHEET> ) {
      @values = split( /,/, $line );
      my $key = $values[$laneColumnIdx].'_';
      my $trueIndex = $values[$indexColumnIdx];
      
      print NANUQ_SHEET $values[0];
      for ( my $idx = 1 ; $idx < @values ; $idx++ ) {
        print NANUQ_SHEET ',';
        if ( $idx == $indexColumnIdx ) {
          if(length($values[$idx]) > 0) {
            my @indexes = split('-', $values[$idx]);
            my $readIdx=0;
            $trueIndex = "";
            for my $rH_readInfo (@{$rAoH_readsInfo}) {
              if($rH_readInfo->{'isIndexed'} eq 'Y') {
                if($readIdx > 0) {
                  print NANUQ_SHEET '-';
                  $key .= '-';
                  $trueIndex .= '-';
                }
                my $index = substr( $indexes[$readIdx], 0,  $rH_readInfo->{'nbCycles'});
                print NANUQ_SHEET $index;
                $key .= $index;
                $trueIndex .= $index;
                $readIdx++;
              }
            }
          }
          else {
            $key .= 'noIndex';
          }
        }
        else {
          print NANUQ_SHEET $values[$idx];
        }
      }
      my %sampleInfo;
      $sampleInfo{'name'} = $values[$sampleIDColumnIdx];
      $sampleInfo{'lane'} = $values[$laneColumnIdx];
      $sampleInfo{'index'} = $trueIndex;
      if(length($sampleInfo{'index'}) == 0 || $sampleInfo{'index'} eq ',') {
        $sampleInfo{'index'} = "NoIndex";
      }
      
      $retVal{$key} = \%sampleInfo;

    }
    close(NANUQ_SHEET);
    close(SSHEET);

  return \%retVal;
}
