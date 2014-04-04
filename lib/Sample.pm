#!/usr/bin/env perl

=head1 NAME

I<Sample>

=head1 SYNOPSIS

Object used to hold information on a Sample

=head1 DESCRIPTION

Object used to hold information on a Sample


=head1 DEPENDENCY

=cut

package Sample;

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
use File::Basename;
use Text::CSV;
use ReadSet;

# SUB
#-----------------------
sub new {
  my $class = shift;
  my $name = shift;

  unless ($name =~ /^\w[\w.-]*$/) {
    die "Error: sample name \"$name\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!";
  }

  my $self = {
    '_name' => $name,
    '_readSets' => [],
  };

  bless($self, $class);
  return $self;
}

sub getName {
  my ($self) = @_;
  return $self->{'_name'};
}

sub getReadSets {
  my ($self) = @_;
  return $self->{'_readSets'};
}

sub getNbReadSets {
  my ($self) = @_;
  return scalar(@{$self->{'_readSets'}});
}

sub addReadSet {
  my ($self, $rO_readSet) = @_;
  if ($self->getReadSetByName($rO_readSet->getName())) {
    die "Error: readSet name \"" . $rO_readSet->getName() . "\" already exists for sample \"" . $self->getName() . "\"!";
  } else {
    push(@{$self->getReadSets()}, $rO_readSet);
    $rO_readSet->setSample($self);
  }
}

sub getReadSetByName {
  my ($self, $name) = @_;

  my @namedReadSets = grep($_->getName() eq $name, @{$self->getReadSets()});

  # Should find 1 readSet at most since readSet names are unique
  if (@namedReadSets) {
    return $namedReadSets[0];
  } else {
    return undef;
  }
}

sub parseSampleFile {
  my ($self, $sampleFile) = @_;

  my $rH_samples = {};

  open my $in, "<$sampleFile" or die "[Error] cannot open $sampleFile: $!";

  my $csv = Text::CSV->new({
    sep_char   => "\t", # TAB-separated values
    auto_diag  => 1, # Report irregularities immediately
    binary     => 1,  # Allow special character.
  });

  $csv->column_names($csv->getline($in)); # Retrieve header

  while (my $inputRow = $csv->getline_hr($in)) {
    my $sampleName = $inputRow->{"Sample"};
    my $sample;

    if ($rH_samples->{$sampleName}) {    # Sample already exists
      $sample = $rH_samples->{$sampleName};
    } else {    # Create new sample
      $sample = Sample->new($sampleName);
      $rH_samples->{$sampleName} = $sample;
    }

    # Create readSet and add it to sample
    my $readSet = ReadSet->new($inputRow->{"ReadSet"}, $inputRow->{"RunType"});

    # ReadSets file paths are either absolute or relative to the sample file
    # Convert them to absolute paths and check if files exist
    foreach my $format ("BAM", "FASTQ1", "FASTQ2") {
      if ($inputRow->{$format}) {
        if (not(File::Spec->file_name_is_absolute($inputRow->{$format}))) {
          $inputRow->{$format} = File::Spec->rel2abs(dirname($sampleFile)) . "/" . $inputRow->{$format};
        }
        -f $inputRow->{$format} or die "[Error] in parseSampleFile: \"" . $inputRow->{$format} . "\" does not exist or is not a valid plain file!";
      }
    }
    $readSet->setBAM($inputRow->{"BAM"});
    $readSet->setFASTQ1($inputRow->{"FASTQ1"});
    $readSet->setFASTQ2($inputRow->{"FASTQ2"});
    $readSet->setLibrary($inputRow->{"Library"});
    $readSet->setRun($inputRow->{"Run"});
    $readSet->setLane($inputRow->{"Lane"});
    $readSet->setAdaptor1($inputRow->{"Adaptor1"});
    $readSet->setAdaptor2($inputRow->{"Adaptor2"});
    $readSet->setQualityOffset($inputRow->{"QualityOffset"});
    $readSet->setBED($inputRow->{"BED"});
    $sample->addReadSet($readSet);
  }

  # Return reference to array of samples
  return [values(%$rH_samples)];
}

1;
