#!/usr/bin/perl -w

use strict;

die "pass a fasta and a fasta-quality file\n"
  unless @ARGV;

my ( $seq_infile, $qual_infile ) =
  ( scalar @ARGV == 1 ) ? ( $ARGV[0], "$ARGV[0].qual" ) : @ARGV;

my $fastaFile;
my $qualFile;
## Create input objects for both a seq (fasta) and qual file

open FASTA, '<', $seq_infile  or die "Couldn't open fasta file.\n";
open QUAL,  '<', $qual_infile or die "Couldn't open quality file.\n";

writeFastq();
close(FASTA);
close(QUAL);

sub writeFastq {
  my $fastaHeader = undef;
  my $qualHeader  = undef;
  while ( my $fastaLine = <FASTA> ) {
    chomp($fastaLine);

    if ( $fastaLine =~ m/^>/ ) {
      $fastaHeader = substr( $fastaLine, 1 );
      last;
    }
  }

  while ( my $qualLine = <QUAL> ) {
    chomp($qualLine);

    if ( $qualLine =~ m/^>/ ) {
      $qualHeader = substr( $qualLine, 1 );
      last;
    }
  }

  if ( $fastaHeader ne $qualHeader ) {
    die "Headers don't match between qual and fasta: $fastaHeader'n";
  }

  my $entryLength = 0;
  while ( my $fastaLine = <FASTA> ) {
    chomp($fastaLine);

    if ( $fastaLine =~ m/^>/ && $entryLength > 0 ) {
      print STDOUT "\n";
      $qualHeader = writeQual( $fastaHeader, $entryLength );
      print STDOUT "\n";
      $fastaHeader = substr( $fastaLine, 1 );
      $entryLength = 0;
      if ( $fastaHeader ne $qualHeader ) {
        die "Headers don't match between qual and fasta: $fastaHeader'n";
      }
    }
    elsif ( defined($fastaHeader) && length( trim($fastaLine) ) != 0 ) {
      if ( $entryLength == 0 ) {
        print STDOUT "@" . $fastaHeader . "\n";
      }

      $entryLength += length($fastaLine);
      print STDOUT $fastaLine;
    }
  }
  print STDOUT "\n";
  writeQual( $fastaHeader, $entryLength );
  print STDOUT "\n";
}

sub writeQual {
  my $fastaHeader = shift;
  my $entryLength = shift;

  my $qualEntryLength = 0;
  print STDOUT "+\n";
  while ( my $qualLine = <QUAL> ) {
    chomp($qualLine);

    if ( $qualLine =~ m/^>/ && $qualEntryLength > 0 ) {
      if ( $qualEntryLength != $entryLength ) {
        die "Entry '"
          . $fastaHeader
          . "' has different fasta vs qual length: "
          . $entryLength . " != "
          . $qualEntryLength . "\n";
      }
      return substr( $qualLine, 1 );
    }
    elsif ( length( trim($qualLine) ) != 0 ) {
      my @values = split( / /, $qualLine );
      foreach my $value (@values) {
        if(length(trim($value)) == 0) {
          next;
        }
        my $ascii = int($value);
        print STDOUT chr( $ascii + 33 );
        $qualEntryLength++;
      }
    }
  }
}

## Trim function.
sub trim {
  my $string = shift;
  $string =~ s/^\s+|\s+$//g;
  return $string;
}
