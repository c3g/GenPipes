#!/usr/bin/perl
use strict;

my $pair1 = $ARGV[0];
my $pair2 = $ARGV[1];
my $out_pair1= $ARGV[2];
my $out_pair2= $ARGV[3];
my $out_single= $ARGV[4];
my $trimLimit = $ARGV[5];

open (PAIR1, "$pair1") || die "Cannot Open pair 1 File";
open (PAIR2, "$pair2") || die "Cannot Open pair 2 File";
open (OUT_PAIR1, '>', "$out_pair1") || die "Cannot Open pair 1 File";
open (OUT_PAIR2, '>', "$out_pair2") || die "Cannot Open pair 2 File";
open (OUT_SINGLE,'>',  "$out_single") || die "Cannot Open single File";

while(my $head1 = <PAIR1>) {
  my $seq1 = <PAIR1>;
  my $qhead1 = <PAIR1>;
  my $qual1 = <PAIR1>;
  my $head2 = <PAIR2>;
  my $seq2 = <PAIR2>;
  my $qhead2 = <PAIR2>;
  my $qual2 = <PAIR2>;

	# We test for +1 because we didn't chomp the input.
	if(length($seq1) < $trimLimit+1 && length($seq2) < $trimLimit+1) {
		next;
	}
	elsif(length($seq1) < $trimLimit+1) {
		print OUT_SINGLE $head2;
		print OUT_SINGLE $seq2;
		print OUT_SINGLE $qhead2;
		print OUT_SINGLE $qual2;
	}
	elsif(length($seq2) < $trimLimit+1) {
		print OUT_SINGLE $head1;
		print OUT_SINGLE $seq1;
		print OUT_SINGLE $qhead1;
		print OUT_SINGLE $qual1;
	}
	else {
		print OUT_PAIR1 $head1;
		print OUT_PAIR1 $seq1;
		print OUT_PAIR1 $qhead1;
		print OUT_PAIR1 $qual1;
		print OUT_PAIR2 $head2;
		print OUT_PAIR2 $seq2;
		print OUT_PAIR2 $qhead2;
		print OUT_PAIR2 $qual2;
	}
}
close(PAIR1);
close(PAIR2);
close(OUT_PAIR1);
close(OUT_PAIR2);
close(OUT_SINGLE);
