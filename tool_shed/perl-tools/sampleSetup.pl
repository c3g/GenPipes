#!/usr/bin/perl

use File::Path qw(mkpath);
use Cwd;
use Text::CSV;
use Getopt::Long;

my $version = "1.0";

use strict;

&main();

sub printUsage {
    print "Usage sampleSetup.pl --nanuqAuthFile \$HOME/.nanuqAuth.txt --usesheet project.nanuq.csv --tech hiseq\n";
    print "\t--nanuqAuthFile <FILE>  Path to nanuq authentication file\n";
    print "\t--usesheet      <FILE>  Use an already existing sample sheet\n";
    print "\t--nolinks               Don't create raw_reads directory or symlinks\n";
    print "\t--projectId     <INT>   Nanuq project id from which to get the sample sheet\n";
    print "\t--help                  This help\n";
    exit(0);
}

sub main {

  my $techName;
  my $projectId;
  my $sampleSheet;
  my $nanuqAuthFile;
  my $noLinks;
  my $help;
  my $result = GetOptions(
    "tech=s"           => \$techName,
    "projectId=s"      => \$projectId,
    "usesheet=s"       => \$sampleSheet,
    "nanuqAuthFilei=s" => \$nanuqAuthFile,
    "nolinks!"         => \$noLinks,
    "help!"         => \$help,
  );

  if($help) {
    printUsage();
  }

  my $errMsg = "";
  if(!defined($nanuqAuthFile) || !-e $nanuqAuthFile) {
    $errMsg .= "Missing nanuqAuthFile\n";
  }
  if(defined($projectId) && defined($sampleSheet)) {
    $errMsg .= "You can't set both projectId and useSheet\n";
  }
  if((!defined($projectId) || length($projectId) == 0) && (!defined($sampleSheet) || length($sampleSheet) == 0)) {
    $errMsg .= "Missing projectId or useSheet\n";
  }
  if(!defined($techName) || length($techName) == 0) {
    $errMsg .= "Missing tech\n";
  }
  if(length($errMsg)) {
    warn $errMsg;
    printUsage();
  }

  my $isMiseq = 0;
  if(lc($techName) eq 'hiseq') {
    $techName = "HiSeq";
  }
  elsif(lc($techName) eq 'miseq') {
    $techName = "MiSeq";
    $isMiseq = 1;
  }

  my $projectFile = 'project.nanuq.csv';
  if(defined($projectId)) {
    getSheet($projectFile, $techName, $projectId, $nanuqAuthFile);
  }
  else {
    $projectFile = $sampleSheet;
  }

  my $rA_SampleInfos = parseSheet($projectFile, $isMiseq);
  if(!$noLinks) {
    createLinks($rA_SampleInfos);
  }
}

sub createLinks {
  my $rA_SampleInfos = shift;

  for my $rH_Sample (@$rA_SampleInfos) {
    my $dirCur = getcwd;
    my @Foldervalues = split('/', $dirCur);
    if ( $Foldervalues[-1] eq "raw_reads" ) {
       my $directory = $rH_Sample->{'name'}."/run".$rH_Sample->{'runId'}."_".$rH_Sample->{'lane'};
    } else {
      my $directory = 'raw_reads/'.$rH_Sample->{'name'}."/run".$rH_Sample->{'runId'}."_".$rH_Sample->{'lane'};
    }
    mkpath($directory);

    my $runType = $rH_Sample->{'runType'};
    my $file1;
    my $file2;
    if($runType eq "SINGLE_END") {
      $file1 = $directory.'/'.$rH_Sample->{'name'}.'.'.$rH_Sample->{'qualOffset'}.".single.fastq.gz";  
    }
    elsif($runType eq "PAIRED_END") {
      $file1 = $directory.'/'.$rH_Sample->{'name'}.'.'.$rH_Sample->{'qualOffset'}.".pair1.fastq.gz";
      $file2 = $directory.'/'.$rH_Sample->{'name'}.'.'.$rH_Sample->{'qualOffset'}.".pair2.fastq.gz";
    }  

    if($runType eq "PAIRED_END") {
      if( ! -l $file1) {
        if(symlink($rH_Sample->{'filename1'}, $file1) or die "Can't symlink ".$rH_Sample->{'filename1'}." -> $file1\n") {
          print "Created link $file1\n";
        }
      }

      if( ! -l $file2) {
        symlink($rH_Sample->{'filename2'}, $file2) or die "Can't symlink ".$rH_Sample->{'filename2'}." -> $file2\n";
      }
    }
    elsif($runType eq "SINGLE_END") {
      if( ! -l $file1) {
        if(symlink($rH_Sample->{'filename1'}, $file1) or die "Can't symlink ".$rH_Sample->{'filename1'}." -> $file1\n") {
          print $rH_Sample->{'filename1'} . "\n";
          print "Created link $file1\n";
        }
      }
    } 
  }
}

sub getSheet {
  my $projectFile = shift;
  my $tech = shift;
  my $projectId = shift;
  my $nanuqAuthFile = shift;

  my $command = 'wget --no-cookies --post-file '.$nanuqAuthFile.' https://genomequebec.mcgill.ca/nanuqMPS/csv/technology/'.$tech.'/project/'.$projectId.'/filename/'.$projectFile."\n";
  print '#'.$command;
  system($command);
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

sub parseSheet {
  my $fileName = shift;
  my  $isMiseq = shift;

  my @retVal;
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";
  my $line = <SAMPLE_SHEET>;
  my $nameIdx=-1;
  my $libraryBarcodeIdx=-1;
  my $runIdIdx=-1;
  my $laneIdx=-1;
  my $runTypeIdx=-1;
  my $statusIdx=-1;
  my $readSetIdIdx=-1;
  my $filePrefixIdx=-1;

  my $csv = Text::CSV->new();
  $csv->parse($line);
  my @headers = $csv->fields();
  for(my $idx=0; $idx < @headers; $idx++) {
    
    $headers[$idx] =~ s/"//g;
    if($headers[$idx] eq "Name") {
      $nameIdx=$idx;
    }
    elsif($headers[$idx] eq "Library Barcode") {
      $libraryBarcodeIdx=$idx;
    }
    elsif($headers[$idx] eq "Run") {
      $runIdIdx=$idx;
    }
    elsif($headers[$idx] eq "Region") {
      $laneIdx=$idx;
    }
    elsif($headers[$idx] eq "Run Type") {
      $runTypeIdx=$idx;
    }
    elsif($headers[$idx] eq "Status") {
      $statusIdx=$idx;
    }
    elsif($headers[$idx] eq "Read Set Id") {
      $readSetIdIdx=$idx;
    }
    elsif($headers[$idx] eq "Filename Prefix") {
      $filePrefixIdx=$idx;
    }
  }

  my $sampleSheetErrors="";
  if($nameIdx==-1) {
    $sampleSheetErrors.="Missing Sample Name\n";
  }
  if($libraryBarcodeIdx==-1) {
    $sampleSheetErrors.="Missing Library Barcode\n";
  }
  if($runIdIdx==-1) {
    $sampleSheetErrors.="Missing Run ID\n";
  }
  if($laneIdx==-1) {
    $sampleSheetErrors.="Missing Lane\n";
  }
  if($runTypeIdx==-1) {
    $sampleSheetErrors.="Missing Run Type\n";
  }
  if($statusIdx==-1) {
    $sampleSheetErrors.="Missing Status\n";
  }
  if($readSetIdIdx==-1) {
    $sampleSheetErrors.="Read Set Id\n";
  }
  if($filePrefixIdx==-1) {
    $sampleSheetErrors.="Filename Prefix\n";
  }
  if(length($sampleSheetErrors) > 0) {
    die $sampleSheetErrors;
  }

  while($line = <SAMPLE_SHEET>) {
    $csv->parse($line);         
    my @values = $csv->fields();
    if($values[$statusIdx] =~ /invalid/) {
      warn "Invalid: $values[$nameIdx] $values[$runIdIdx] $values[$laneIdx]\n";
      next;
    }

    my %sampleInfo;
    $sampleInfo{'name'} = $values[$nameIdx];
    $sampleInfo{'libraryBarcode'} = $values[$libraryBarcodeIdx];
    $sampleInfo{'runId'} = $values[$runIdIdx];
    $sampleInfo{'lane'} = $values[$laneIdx];
    $sampleInfo{'runType'} = $values[$runTypeIdx];
    $sampleInfo{'readSetId'} = $values[$readSetIdIdx];
    $sampleInfo{'filePrefix'} = $values[$filePrefixIdx];

    my $rootDir;
    my $isHiSeq = 0;
    if($isMiseq == 0) {
      $rootDir = "/lb/robot/hiSeqSequencer/hiSeqRuns/";
      $isHiSeq = 1;
    }
    elsif($isMiseq == 1) {
      $rootDir = "/lb/robot/miSeqSequencer/miSeqRuns/";
    }
    else {
      die "Unknown prefix technology type: ".$sampleInfo{'filePrefix'}."\n";
    }
    opendir(ROOT_DIR, $rootDir) or die "Couldn't open directory ".$rootDir."\n";
    my @rootFiles;
    if($isHiSeq == 1) {
      @rootFiles =  grep { /.*[0-9]+_[^_]+_[^_]+_$sampleInfo{'runId'}/ } readdir(ROOT_DIR);
      if(@rootFiles == 0) {
        opendir(ROOT_DIR_2013, $rootDir .'/2013') or die "Couldn't open directory ".$rootDir."/2013\n";
        rootFiles =  grep { /.*[0-9]+_[^_]+_[^_]+_$sampleInfo{'runId'}/ } readdir(ROOT_DIR_2013);
      }
    }
    else {
      @rootFiles =  grep { /.*[0-9]+_$sampleInfo{'runId'}/ } readdir(ROOT_DIR);
      if(@rootFiles == 0) {
        opendir(ROOT_DIR_2013, $rootDir .'/2013') or die "Couldn't open directory ".$rootDir."/2013\n";
        rootFiles =  grep { /.*[0-9]+_$sampleInfo{'runId'}/ } readdir(ROOT_DIR_2013);
      }
    }

    if(@rootFiles == 0) {
      die "Run not found: ".$sampleInfo{'runId'}."\n";
    }
    elsif(@rootFiles > 1) {
      die "Many runs found: ".$sampleInfo{'runId'}."\n";
    }

    my $runPath  = $rootDir.'/'.$rootFiles[0];
    my $fastqDir = `echo $runPath/se*`;
    chomp($fastqDir);
    if($fastqDir =~ /\*/){
      $fastqDir = `echo $runPath/Data/In*/B*/G*`;
      chomp($fastqDir);
      if($fastqDir =~ /\*/){
        die "Couldn't find fastq directory: $fastqDir\n";
      }
      $sampleInfo{'qualOffset'} = "64";
      my $toTest1 = $fastqDir.'/s_'.$sampleInfo{'lane'}.'_1_*'.$sampleInfo{'name'}.'*.txt.gz';
      my $toTest2 = $fastqDir.'/s_'.$sampleInfo{'lane'}.'_2_*'.$sampleInfo{'name'}.'*.txt.gz';
      $sampleInfo{'filename1'} = `echo $toTest1`;
      $sampleInfo{'filename2'} = `echo $toTest2`;
    }
    else {
      $sampleInfo{'qualOffset'} = "33";
     
      my $toTest1;
      my $toTest2;
      my $typeRun = $values[$runTypeIdx];
      if($typeRun eq "PAIRED_END") {
        $toTest1 = $fastqDir.'/'.$sampleInfo{'filePrefix'}.'_R1.fastq.gz';
        $toTest2 = $fastqDir.'/'.$sampleInfo{'filePrefix'}.'_R2.fastq.gz';
        $sampleInfo{'filename1'} = `echo $toTest1`;
        $sampleInfo{'filename2'} = `echo $toTest2`;
      }
      elsif($typeRun eq "SINGLE_END") {
        $toTest1 = $fastqDir.'/'.$sampleInfo{'filePrefix'}.'_R1.fastq.gz';
        $sampleInfo{'filename1'} = `echo $toTest1`;
      } 
    }
    
    chomp($sampleInfo{'filename1'});
    if($values[$runTypeIdx] eq "PAIRED_END") {
      chomp($sampleInfo{'filename2'});
    }
    push(@retVal, \%sampleInfo);
  }

  return \@retVal;
}

