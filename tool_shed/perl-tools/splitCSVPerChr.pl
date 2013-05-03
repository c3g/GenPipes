#!/usr/bin/perl

use strict;


&main();

sub main {
	my $fileToSplit = $ARGV[0];

	open(FILE, $fileToSplit) or die "Can't open $fileToSplit\n";
  my ($newFileName) = $fileToSplit =~ /(.+)\.[^.]+$/g;
	my $header = <FILE>;
	my $currentChromosome = "-1";
	my $newFH = undef;
  my $fileIdx=1;
  my $lineCount=0;
	while(my $line = <FILE>) {
		my @values = split(/\t/, $line);
    $lineCount++;

		if($currentChromosome != $values[0]) {
			if(defined $newFH) {
				close($newFH);
			}
      $currentChromosome = $values[0];
      $fileIdx=1;
			open($newFH, ">".$newFileName.'.'.$currentChromosome.'.'.$fileIdx.'.tsv') or die "Can't write to ".$fileToSplit.'.'.$currentChromosome."\n";
      $lineCount=0;
			print { $newFH } $header;
		}
    elsif($lineCount > 100000) {
      close($newFH);
      $fileIdx++;
      open($newFH, ">".$newFileName.'.'.$currentChromosome.'.'.$fileIdx.'.tsv') or die "Can't write to ".$fileToSplit.'.'.$currentChromosome."\n";
      $lineCount=0;
      print { $newFH } $header;
    }
		print { $newFH } $line;
	}
	if(defined $newFH) {
    close($newFH);
  }
	close(FILE);
}
