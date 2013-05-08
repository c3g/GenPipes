#!/usr/bin/perl
use Time::HiRes qw( time );

use strict;
use warnings;
use Getopt::Std;

use vars qw($VERSION);
$VERSION = '1.0';

&main;

sub main {

  my %opts = ( r => undef, o => 'output.out' );
  getopts( 'r:o:', \%opts );

  if ( !defined( $opts{'r'} ) ) {

    die(
    qq/
Usage:   readsPerOasesTranscript.pl [options] -r <Velvet\/Oases run directory> -o <output>
Version: $VERSION
Options: -r        Velvet\/Oases run directory (Directory with LastGraph, Sequences, contig-ordering, Log files).
         -o File   Output file [$opts{o}]
\n/
    )
  }

  my $lastGraph = $opts{'r'}."\/LastGraph";
  my $sequences = $opts{'r'}."\/Sequences";
  my $contigOrdering = $opts{'r'}."\/contig-ordering.txt";

  my $rA_readIdRanges = getReadIdRanges( $sequences );

  #  for my $ranges (@$rA_readIdRanges) {
  #    print "Ranges: $ranges\n";
  #  }

  my $rHoA_contigs = parseContigsFromGraph( $lastGraph, $rA_readIdRanges );
  my $rHoA_transcripts = parseContigOrder( $contigOrdering, $rHoA_contigs );

  for my $transcriptName ( keys %{$rHoA_transcripts} ) {
    print $transcriptName;
    for my $value ( @{ $rHoA_transcripts->{$transcriptName} } ) {
      print "\t";
      print $value;
    }
    print "\n";
  }
}

sub getReadIdRanges {
  my $sequencesFile = shift;

  my @readIdRanges;
  my $prevReadId;
  my $prevHeader;
  my $prevGroup = -1;
  my $header;
  my $readId;
  my $readGroup;
  open( SEQ_FILE, "$sequencesFile" ) or die "Can't open $sequencesFile\n";
  while (<SEQ_FILE>) {
    chomp;
    my $line = $_;
    if ( $line =~ /^>/ ) {
      ( $header, $readId, $readGroup ) = split(/\t/);
      if ( $prevGroup != $readGroup ) {
        if ( $prevGroup != -1 ) {
          print STDERR "$prevHeader\t$prevReadId\n";
          push( @readIdRanges, $prevReadId );
        }
        $prevGroup = $readGroup;
      }
      $prevReadId = $readId;
      $prevHeader = $header;
    }
  }
  # Add the last one
  print STDERR "$header\t$readId\n";
  push( @readIdRanges, $readId );
  close(SEQ_FILE);

  return \@readIdRanges;
}

sub parseContigOrder {
  my $contigOrderFile = shift;
  my $rHoA_contigs    = shift;

  print STDERR "Started to parse Contig Order file\n";
  my $start = time();

  my %transcriptInfo;
  open( CONTIG_ORDER, $contigOrderFile )
    or die "Couldn't open contig-order file: $contigOrderFile\n";
  while ( my $line = <CONTIG_ORDER> ) {
    chomp($line);
    if ( $line =~ /^>(.*Transcript.*)/ ) {
      my $transcriptName = $1;

      $line = <CONTIG_ORDER>;
      chomp($line);

      #$contig_id:$cumulative_length-($distance_to_next_contig)->
      my @contigDetailList = split( /->/, $line );

      my @transcriptReadCount;
      for my $contigDetail (@contigDetailList) {
        my ($contigId) = $contigDetail =~ /^([^:]+):/;

        my $foundIt = 0;
        if ( !defined $rHoA_contigs->{$contigId} ) {
          $contigId = $contigId * -1;
          if ( !defined $rHoA_contigs->{$contigId} ) {
            warn( "Contig [" . ( $contigId ) . "] wasn't found\n" );
          }
          else {
            $foundIt = 1;
          }
        }
        else {
          $foundIt = 1;
        }

        if ($foundIt) {
          for ( my $idx = 0 ; $idx < @{ $rHoA_contigs->{$contigId} } ; $idx++ )
          {
            $transcriptReadCount[$idx] += $rHoA_contigs->{$contigId}->[$idx];
          }
        }
      }
      $transcriptInfo{$transcriptName} = \@transcriptReadCount;
    }
  }
  close(CONTIG_ORDER);

  my $end = time();
  printf STDERR ( "[%.2f] Parsed Contig Order file\n", $end - $start );
  return \%transcriptInfo;
}

sub parseContigsFromGraph {
  my $graphFile       = shift;
  my $rA_readIdRanges = shift;

  my @nbReadsUsedPerRange;

  print STDERR "Started to parse Graph file\n";
  my $start = time();
  my %contigs;
  my $graphFilesize = -s $graphFile;
  open( GRAPH, $graphFile ) or die "Couldn't open graph file: $graphFile\n";
  my $readingShortRead       = 0;
  my $readingLongReadContigs = 0;
  print STDERR "\n";

  while ( my $line = <GRAPH> ) {
    chomp $line;
    my $currentPos = tell(GRAPH);
#    if ( ( $. % 100000 ) == 0 ) {
#      print STDERR $currentPos . "/" . $graphFilesize . "\r";
#    }

    # Short read line
    if ( $line =~ /^NR\t([^\t]+)\t([0-9]+)/ ) {
      $readingShortRead       = $1;
      $readingLongReadContigs = 0;
      if ( !defined $contigs{$1} ) {
        $contigs{$1} = [];
        for ( my $idx = 0 ; $idx < @$rA_readIdRanges ; $idx++ ) {
          $contigs{$1}->[$idx] = 0;
        }
      }
    }
    # Long read line
    elsif ( $line =~ /^SEQ\t-?(.*)/ ) {
      my $readId = 0 + $1;
      my $idx = 0;
      my $foundSource = 0;
      for ( $idx = 0 ; $idx < @$rA_readIdRanges ; $idx++ ) {
        if ( $readId < $rA_readIdRanges->[$idx] ) {
          $foundSource = 1;
          last;
        }
      }

      if($foundSource == 0) {
        warn( "Found offseted read: " . $readId . "\n" );
      }
      elsif ($idx != (@{$rA_readIdRanges}-1)) {
        warn( "Long read not the last found in Sequences: " . $readId . "\n" );
      }
      if(! defined($nbReadsUsedPerRange[$idx])) {
				$nbReadsUsedPerRange[$idx] = {};
			}
      $nbReadsUsedPerRange[$idx]->{$readId} = 1;
      
      $readingShortRead       = 0;
      $readingLongReadContigs = 1;
    }
    elsif ( $readingLongReadContigs || $readingShortRead ) {
      if ( $line =~ /^(-[0-9]+)/ || $line =~ /^([0-9]+)/ ) {
        if ($readingLongReadContigs) {
          my $contigId = $1;
					print STDERR "Reading long contig: ".$contigId."\n";

          if ( !defined $contigs{$contigId} ) {
            $contigs{$contigId} = [];
            for ( my $idx = 0 ; $idx < @$rA_readIdRanges ; $idx++ ) {
              $contigs{$contigId}->[$idx] = 0;
            }
          }
          my $value = $contigs{$contigId}->[ @{$rA_readIdRanges} - 1 ];
          $value++;
          $contigs{$contigId}->[ @{$rA_readIdRanges} - 1 ] = $value;
					print STDERR "Reading long contig: ".$contigId." value ".$value."\n";
        }
        elsif ($readingShortRead) {
          my $readId = 0 + $1;
          $readId = 0 + $readId;

          #$usedReads{$readId} = 1;

          my $readAdded = 0;
          my $idx = 0;
          for ( $idx = 0 ; $idx < @$rA_readIdRanges ; $idx++ ) {
            if ( $readId < $rA_readIdRanges->[$idx] ) {
              my $value = $contigs{$readingShortRead}->[$idx];
              $value++;
              $contigs{$readingShortRead}->[$idx] = $value;
              $readAdded = 1;
              last;
            }
          }

      		if(! defined($nbReadsUsedPerRange[$idx])) {
						$nbReadsUsedPerRange[$idx] = {};
					}
		      $nbReadsUsedPerRange[$idx]->{$readId} = 1;

          if($readAdded == 0) {
            warn( "Found offseted read: " . $readId . "\n" );
          }
        }
      }
      else {
        $readingShortRead = 0;
        $readingLongReadContigs = 0;
      }
    }
  }
  close(GRAPH);

  my $end = time();
  printf STDERR ( "[%.2f] Parsed graph file\n", $end - $start );

  for my $rH_readIds (@nbReadsUsedPerRange) {
	  print "\t".keys(%{$rH_readIds})
  }
  print "\n";

  #for my $rr (keys %usedReads) {
  #print STDERR $rr."\n";
  #}
  #print STDERR "done\n";
  return \%contigs;
}

