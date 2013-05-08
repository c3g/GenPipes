#!/usr/bin/perl

use Getopt::Long;
use strict;

if ( $#ARGV < 0 ) {
  print STDERR <<EOF;
Options:
--qual-threshold <num>     quality threshold for trimming, default 20
--length-threshold <num>    length threshold for trimming, default 20
--qual-type <num>           0=sanger qualities, 1=illumina qualities pipeline>=1.3, 2=illumina qualities pipeline<1.3.  Default 0.
--pair1 <paired end input filename>      fastq, paired end file. Must have same number of records as pair2. Required.
--pair2 <paired end input filename>      fastq, paired end file. Must have same number of records as pair1. Required.
--outpair1 <paired end output file>      Required.
--outpair2 <paired end output file>      Required.
--single <single end output file>        Required.
EOF
  exit(1);
}

my $qual_threshold   = 20;
my $length_threshold = 20;
my $qual_type        = 0;
my ( $pair1, $pair2, $outpair1, $outpair2, $single );

GetOptions(
  'qual-threshold=i'   => \$qual_threshold,
  'length-threshold=i' => \$length_threshold,
  'qual-type=i'        => \$qual_type,
  'pair1=s'            => \$pair1,
  'pair2=s'            => \$pair2,
  'outpair1=s'         => \$outpair1,
  'outpair2=s'         => \$outpair2,
  'single=s'           => \$single
);

if ( !defined $pair1
  || !defined $pair2
  || !defined $outpair1
  || !defined $outpair2
  || !defined $single )
{
  print STDERR
    "Error: You must define pair1, pair2, outpair1, outpair2, and single\n";
  exit(1);
}

sub get_quality_num {
  my ( $qual_char, $qual_type ) = @_;

  if ( $qual_type == 0 ) {
    return ( ord($qual_char) - 33 );
  }

  elsif ( $qual_type == 1 ) {
    return ( ord($qual_char) - 64 );
  }

  elsif ( $qual_type == 2 ) {
    return (
      10 * log( 1 + 10**( ( ord($qual_char) - 64 ) / 10.0 ) ) / log(10) );
  }
}

open( PAIR1, '-|', 'zcat', "$pair1" ) or die "Cannot access file $pair1.\n";
open( PAIR2, '-|', 'zcat', "$pair2" ) or die "Cannot access file $pair2.\n";
open( OUTPAIR1, ">$outpair1" );
open( OUTPAIR2, ">$outpair2" );
open( SINGLE,   ">$single" );

my (
  $count,      $p1_header1, $p1_seq,     $p1_header2, $p1_qual,
  $p2_header1, $p2_seq,     $p2_header2, $p2_qual
);
my ( $pos, $p1_flag, $qual_char, $qualnum, $p1_cut, $p2_flag, $p2_cut );
my ( $paired_cnt, $single_cnt, $discard_cnt ) = ( 0, 0, 0 );
my (
  $i,            $window_total, $window_size,
  $window_start, @p1_qual_data, @p2_qual_data
);
my ( $trimmed_p1_seq, $trimmed_p1_qual, $trimmed_p2_seq, $trimmed_p2_qual );

$count = 0;
while ( $p1_header1 = <PAIR1> ) {
  $p1_seq     = <PAIR1>;
  $p1_header2 = <PAIR1>;
  $p1_qual    = <PAIR1>;
  chomp $p1_seq;
  chomp $p1_qual;

  $p2_header1 = <PAIR2>;
  $p2_seq     = <PAIR2>;
  $p2_header2 = <PAIR2>;
  $p2_qual    = <PAIR2>;
  chomp $p2_seq;
  chomp $p2_qual;

  if (
    substr( $p1_header1, 0, length($p1_header1) - 2 ) ne
    substr( $p2_header1, 0, length($p2_header1) - 2 ) )
  {
    die("Headers don't match. 1: $p1_header1, 2: $p2_header1\n");
  }

  $pos     = 0;
  $p1_flag = 0;
  my $rev_qual = reverse($p1_qual);
  foreach $qual_char ( split( //, $rev_qual ) ) {
    $qualnum = get_quality_num( $qual_char, $qual_type );

    if ( $qualnum >= $qual_threshold ) {
      $p1_cut = length($p1_qual) - $pos;
      if ( $p1_cut >= $length_threshold ) {
        $p1_flag = 1;
      }

      last;
    }

    $pos++;
  }

  $pos      = 0;
  $p2_flag  = 0;
  $rev_qual = reverse($p2_qual);
  foreach $qual_char ( split( //, $rev_qual ) ) {
    $qualnum = get_quality_num( $qual_char, $qual_type );

    if ( $qualnum >= $qual_threshold ) {
      $p2_cut = length($p2_qual) - $pos;
      if ( $p2_cut >= $length_threshold ) {
        $p2_flag = 1;
      }
      last;
    }

    $pos++;
  }

  if ( $p1_flag == 1 && $p2_flag == 1 ) {
    $trimmed_p1_seq  = substr( $p1_seq,  0, $p1_cut );
    $trimmed_p1_qual = substr( $p1_qual, 0, $p1_cut );
    print OUTPAIR1 "$p1_header1$trimmed_p1_seq\n$p1_header2$trimmed_p1_qual\n";

    $trimmed_p2_seq  = substr( $p2_seq,  0, $p2_cut );
    $trimmed_p2_qual = substr( $p2_qual, 0, $p2_cut );
    print OUTPAIR2 "$p2_header1$trimmed_p2_seq\n$p2_header2$trimmed_p2_qual\n";

    $paired_cnt++;
  }

  elsif ( $p1_flag == 1 && $p2_flag == 0 ) {
    $trimmed_p1_seq  = substr( $p1_seq,  0, $p1_cut );
    $trimmed_p1_qual = substr( $p1_qual, 0, $p1_cut );
    print SINGLE "$p1_header1$trimmed_p1_seq\n$p1_header2$trimmed_p1_qual\n";

    $single_cnt++;
  }

  elsif ( $p1_flag == 0 && $p2_flag == 1 ) {
    $trimmed_p2_seq  = substr( $p2_seq,  0, $p2_cut );
    $trimmed_p2_qual = substr( $p2_qual, 0, $p2_cut );
    print SINGLE "$p2_header1$trimmed_p2_seq\n$p2_header2$trimmed_p2_qual\n";

    $single_cnt++;
  }

  else { $discard_cnt++; }

  $count++;
  if ( $count % 10000 == 0 ) {
    print STDERR
"Records: $count, Paired: $paired_cnt, Single: $single_cnt, Discarded: $discard_cnt                                                         \r";
  }
}

print STDERR
"Records: $count, Paired: $paired_cnt, Single: $single_cnt, Discarded: $discard_cnt                                                         \n";

close(PAIR1);
close(PAIR2);
close(OUTPAIR1);
close(OUTPAIR2);
close(SINGLE);

