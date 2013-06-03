#!/usr/bin/perl -w

use strict;

die "pass a fasta and a fasta-quality file\n"
  unless @ARGV;

my ( $fastqFile, $fastaFile, $qualFile ) = @ARGV;

open FASTQ, '<', $fastqFile  or die "Couldn't open fastq file.\n";
open FASTA, '>', $fastaFile  or die "Couldn't write to fasta file.\n";
open QUAL,  '>', $qualFile or die "Couldn't write to quality file.\n";

writeFastaQual();
close(FASTA);
close(QUAL);
close(FASTQ);

sub writeFastaQual {
  my $fastqHeader = undef;
  
  # State values
  # 0- Find header
  # 1- Read fasta / Find Qual
  # 2- Read Qual
  my $state = 0;
  my $fastaLength = 0;
  my $qualityLength = 0;
  while ( my $fastqLine = <FASTQ> ) {
    chomp($fastqLine);

    if ( $state == 0 && $fastqLine =~ m/^@/ ) {
      if(defined($fastqHeader)) {
        print FASTA "\n";
      }
      
      # Read header
      $fastqHeader = substr( $fastqLine, 1 );
      print FASTA ">".$fastqHeader."\n";
      $fastaLength = 0;
      $state = 1
    }
    elsif ($state == 1 && $fastqLine =~ m/^\+/ ) {
      # Qual starts
      print FASTA "\n";
      print QUAL ">".$fastqHeader."\n";
      $qualityLength = 0;
      $state = 2;
    }
    elsif ($state == 1) {
      print FASTA $fastqLine;
      $fastaLength += length($fastqLine);
    }
    elsif ($state == 2) {
      my @quals = split(//, $fastqLine);
      foreach my $fastqQual (@quals) {
        $qualityLength++;
        if($qualityLength != 0) {
          print QUAL " ";
        }
        my $qual = ord($fastqQual) - 33;
        if($qual < 0) {
          die "Wrong qual format, not phred+33\n";
        }
        print QUAL $qual;
      }
      if($qualityLength > $fastaLength) {
        die "More qualities than fasta data was found.\n";
      }
      elsif($qualityLength == $fastaLength) {
        $state = 0;
        print QUAL "\n";
      } 
    }
  }
}
