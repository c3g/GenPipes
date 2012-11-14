#!/usr/bin/perl

=head1 NAME

I<SampleSheet>

=head1 SYNOPSIS

SampleSheet->parseSampleSheet(sampleSheet_file_name)

=head1 DESCRIPTION

B<SampleSheet> is a library that parses a Nanuq generated
sample sheet and populates an array with he parsed values.
Each row from the sheet is a new hash in the array.

Input = /path/samplesheet_file_name

Output = %hash 


=head1 AUTHOR

B<Louis Letourneau> - I<louis.letourneau@mail.mcgill.ca>
B<Maxime Caron> - I<max.caron@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<File::Basename> path parsing

B<Cwd> path parsing

=cut

package SampleSheet;

# Strict Pragmas
#---------------------
use strict;
use warnings;
#---------------------

# Dependencies
#--------------------
use File::Basename;
use Cwd 'abs_path';
#--------------------


# SUB
#--------------------

sub parseSampleSheetAsHash {
  my $fileName = shift;

  my $rA_SampleLaneInfos = parseSampleSheet($fileName);
  my %sampleInfo;
  for my $rH_Sample (@$rA_SampleLaneInfos) {
    if(!defined $sampleInfo{ $rH_Sample->{'name'} }) {
      $sampleInfo{ $rH_Sample->{'name'} } = [];
    }

    push(@{$sampleInfo{ $rH_Sample->{'name'} }}, $rH_Sample);
  }
  return \%sampleInfo;
}

sub parseSampleSheet {
  my $fileName = shift;

  my @retVal;
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";
  my $line = <SAMPLE_SHEET>;
  my @headers = split(/",/, $line);
  my ($nameIdx,$libraryBarcodeIdx,$runIdIdx,$laneIdx,$runTypeIdx,$statusIdx,$qualOffsetIdx) = parseHeaderIndexes(\@headers);

  while($line = <SAMPLE_SHEET>) {
    $line =~ s/"//g;
    my @values = split(/,/, $line);
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
    $sampleInfo{'qualOffset'} = $values[$qualOffsetIdx];

    if($values[$runTypeIdx] eq "PAIRED_END") {
      $sampleInfo{'read1File'} = $sampleInfo{'name'}.'.'.$sampleInfo{'qualOffset'}.".pair1.fastq.gz";
      $sampleInfo{'read2File'} = $sampleInfo{'name'}.'.'.$sampleInfo{'qualOffset'}.".pair2.fastq.gz";
    }
    elsif($values[4] eq "SINGLE_END") {
      $sampleInfo{'read1File'} = $sampleInfo{'name'}.'.'.$sampleInfo{'qualOffset'}.".single.fastq.gz";
    }
    else {
      print "Unrecognized run type $values[$runTypeIdx] \n";
      exit 1;
    }

    push(@retVal, \%sampleInfo);
  }

  return \@retVal;
}

sub parseHeaderIndexes {
	my $rA_headers = shift;
	my $nameIdx=-1;
  my $libraryBarcodeIdx=-1;
  my $runIdIdx=-1;
  my $laneIdx=-1;
  my $runTypeIdx=-1;
  my $statusIdx=-1;
  my $qualOffsetIdx=-1;
	
	for(my $idx=0; $idx < @{$rA_headers}; $idx++) {
		my $header = $rA_headers->[$idx];
    $header =~ s/"//g;
    if($header eq "Name") {
      $nameIdx=$idx;
    }
    elsif($header eq "Library Barcode") {
      $libraryBarcodeIdx=$idx;
    }
    elsif($header eq "Run") {
      $runIdIdx=$idx;
    }
    elsif($header eq "Region") {
      $laneIdx=$idx;
    }
    elsif($header eq "Run Type") {
      $runTypeIdx=$idx;
    }
    elsif($header eq "Status") {
      $statusIdx=$idx;
    }
    elsif($header eq "Quality Offset") {
      $qualOffsetIdx=$idx;
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
  if($qualOffsetIdx==-1) {
      $sampleSheetErrors.="Missing Quality Offset\n";
    }
  
  if(length($sampleSheetErrors) > 0) {
    die $sampleSheetErrors;
  }
  
  return ($nameIdx,$libraryBarcodeIdx,$runIdIdx,$laneIdx,$runTypeIdx,$statusIdx,$qualOffsetIdx);
}
1;
