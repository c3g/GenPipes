#!/usr/bin/perl

use strict;

&main;

sub main {
  my $vcfFile = $ARGV[0];

  open (VCF, $vcfFile) or die "Can't open $vcfFile\n";
  while(my $line = <VCF>) {
    chomp($line);
    my @values = split('\t', $line);
# $values[4] == ALT
    if($values[3] =~ /N{1000,}/) {
      print STDERR $line."\n";
    }
    else {
      print $line."\n";
    }
  }
  close(VCF);
}
