#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

use vars qw($VERSION);
$VERSION = '1.0';

&main();

sub main() {

	my %opts = (1=>undef, 2=>undef, o=>'output.out');
  getopts('1:2:o:', \%opts);

  die(qq/
Usage:   pileupCmp.pl [options] -1 <pileup_file_1> -2 <pileup_file_2>
Version: $VERSION
Options: -1        First pileup file
         -2        Second pileup file
         -o File   Output file [$opts{o}]
\n/) if (!defined($opts{1}) || !defined($opts{2}));

	open(FILE1, "<$opts{1}") or die "Cannot open file $opts{1}\n";
	open(FILE2, "<$opts{2}") or die "Cannot open file $opts{2}\n";

	my %chr2pos;
	while(my $line = <FILE1>) {
		chomp($line);
		my @values = split(/\t/, $line);
		if(!defined($chr2pos{$values[0]})) {
			$chr2pos{$values[0]} = {};
		}
		$chr2pos{$values[0]}->{$values[1]} = $values[3];
	}
	close(FILE1);

	
	open(OUTPUT, ">$opts{o}") or die "Cannot create file $opts{o}\n";
	my $identicalVariations = 0;
	while(my $line = <FILE2>) {
		chomp($line);
		my @values = split(/\t/, $line);

		if(defined($chr2pos{$values[0]})) {
			if(defined($chr2pos{$values[0]}->{$values[1]}) && $chr2pos{$values[0]}->{$values[1]} eq $values[3]) {
				$identicalVariations++;
			}
			else {
				print OUTPUT "$line\n";
			}
		} else {
				print OUTPUT "$line\n";
		}
	}
	close(OUTPUT);
	close(FILE2);
	print STDOUT "Identical snps: $identicalVariations\n";
}
