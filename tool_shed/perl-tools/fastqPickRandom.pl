#!/usr/bin/perl

use Getopt::Long;
use strict;

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

if ($#ARGV < 0) {
    print STDERR <<EOF;
Options:
--threshold <num>                       0..1 Random picker default 0.5
--input1 <input filename>               fastq file. If input2 is given this file Must have same number of records as input2. Required.
--input2 <paired end input filename>    fastq, paired end file. Must have same number of records as input1. Required.
--out1 <paired end output file>         Required.
--out2 <paired end output file>         Output file
--compressed                            Use compressed fastqs
--validationPattern <regexp>            Pair regular expression. Pass nothing to turn off validation
EOF
    exit(1);
}

my $threshold = 0.5;
my ($input1, $input2, $out1, $out2, $compressed, $validationPattern);

GetOptions ('threshold=f' => \$threshold,
            'input1=s' => \$input1,
            'input2=s' => \$input2,
            'out1=s' => \$out1,
            'out2=s' => \$out2,
            'compressed!' => \$compressed,
            'validationPattern=s' => \$validationPattern);

if (!defined $input1 || !defined $out1) {
    print STDERR "Error: You must define input1 and out1\n";
    exit(1);
}
print STDERR "Compressed: $compressed\n";
if($compressed) {
print STDERR "Compressed: Y\n";
  
  open (INPUT1, '-|', 'zcat', "$input1") or die "Cannot access file $input1.\n";
  if(defined($input2)) {
    open (INPUT2, '-|', 'zcat', "$input2") or die "Cannot access file $input2.\n";
  }
}
else {
print STDERR "Compressed: N\n";
  open (INPUT1, "$input1") or die "Cannot access file $input1.\n";
  if(defined($input2)) {
    open (INPUT2, "$input2") or die "Cannot access file $input2.\n";
  }
}
open (OUT1, ">$out1");
if(defined($input2)) {
  open (OUT2, ">$out2");
}

my ($count, $p1_header1, $p1_seq, $p1_header2, $p1_qual, $p2_header1, $p2_seq, $p2_header2, $p2_qual);
my $hit_cnt = 0;

$count = 0;
while ($p1_header1 = <INPUT1>) {
    $p1_seq = <INPUT1>;
    $p1_header2 = <INPUT1>;
    $p1_qual = <INPUT1>;
    chomp $p1_seq;
    chomp $p1_qual;
    if(defined($input2)) {
      $p2_header1 = <INPUT2>;
      $p2_seq = <INPUT2>;
      $p2_header2 = <INPUT2>;
      $p2_qual = <INPUT2>;
      chomp $p2_seq;
      chomp $p2_qual;
    }

    if(!defined($validationPattern)) {
      if($p1_header1 =~ /\/1$/){
				$validationPattern = "(.+)\/\\d\$";
			}
      elsif($p1_header1 =~ /\d:[YN]:\d:[ATCG]*$/){
				$validationPattern = "(.+) \\d:[YN]:\\d:[ATCG]*\$";
			}
      else {
        $validationPattern = "";
      }
		}

    if(length(trim($validationPattern)) > 0) {
      $p1_header1 =~ m/$validationPattern/;
      my $head1 = $1;
      if(defined($input2)) {
        $p2_header1 =~ m/$validationPattern/;
        my $head2 = $1;

        if($head1 ne $head2) {
          die ("Headers don't match. 1: $p1_header1, 2: $p2_header1\n");
        }
      }
    }

    my $dice = rand();
    if($dice <= $threshold) {
      $hit_cnt++;
      print OUT1 "$p1_header1$p1_seq\n$p1_header2$p1_qual\n";
      if(defined($input2)) {
        print OUT2 "$p2_header1$p2_seq\n$p2_header2$p2_qual\n";
      }
    }

    $count++;
    if ($count % 10000 == 0) {
      print STDERR "Records: $count, Hits: $hit_cnt                                                         \r";
    }
}

print STDERR "Records: $count, Hits: $hit_cnt                                                         \n";
print STDERR "Hit Percent: ".($hit_cnt*100/$count)."\n";

close (OUT1);close (INPUT1);
if(defined($input2)) {
  close (OUT2);
  close (INPUT2);
}



