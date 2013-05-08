#!/usr/bin/perl -w

use strict;
use Bio::SearchIO;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::DB::EUtilities;

&main();

sub main {
  my $rAoH_transcripts = blast($ARGV[0]);

  my $length = scalar(@$rAoH_transcripts);
  print "Number of transcripts: ".$length."\n";

  print "Name\tHit Accession\tHit eScore\tHit description\t[Other Hit Accession\tOther Hit eScore\tOther Hit description]\n";
  for my $rH_transcript (@$rAoH_transcripts) {
    print $rH_transcript->{'name'};
    print "\t";
		if(!defined $rH_transcript->{'hitAccession'}) {
			print "No Hits";
		}
		else {
  	  print $rH_transcript->{'hitAccession'};
	    print "\t";
    	print $rH_transcript->{'hitSignificance'};
  	  print "\t";
	    print $rH_transcript->{'hitLength'};
    	print "\t";
  	  print $rH_transcript->{'hitDescription'};
	    for my $rH_otherHit (@{$rH_transcript->{'otherHits'}}) {
    	  print "\t";
  	    print $rH_otherHit->{'hitAccession'};
	      print "\t";
      	print $rH_otherHit->{'hitSignificance'};
    	  print "\t";
  	    print $rH_otherHit->{'hitLength'};
	    }
		}
    print "\n";
  }
}

sub blast {
  my $blastDirectory =shift;
  my @transcripts;

  #opendir DIR, $blastDirectory or die "cannot open dir $blastDirectory: $!";
  #my @files = grep { /\.xml$/ } readdir DIR;
  #closedir DIR;
  my @files = ($blastDirectory);

  my $fileCount = scalar(@files);
  my $filesDone = 0;
  print STDERR "\n";
  for my $blastOutput (@files) {
    my $in = new Bio::SearchIO(
      -format => 'blastxml',
#      -file   => $blastDirectory.$blastOutput
      -file   => $blastOutput
    );

    while ( my $result = $in->next_result ) {
      my @query = split( / /, $result->query_description, );
      my $queryName = substr( $query[0], 0, length( $query[0] ));

      my %transcript;
      $transcript{'name'} = $queryName;

      my @stats  = $result->available_statistics;
      my @params = $result->available_parameters;

      my $first = 1;
      while (my $hit = $result->next_hit) {
        $hit->name =~ /gi\|(\d+)/;
        if ($first) {
          $transcript{'hitGI'}              = $1;
          $transcript{'hitAccession'}       = $hit->accession;
          $transcript{'hitFullDescription'} = $hit->description;
          $transcript{'hitSignificance'}    = $hit->significance;
          $transcript{'hitScore'}           = $hit->raw_score;
          $transcript{'hitLength'}          = $hit->length;

          my @entries = split( /\|/, $transcript{'hitFullDescription'} );
          $transcript{'hitDescription'} = $entries[0];

          $transcript{'otherHits'}      = [];
          $first = 0;
        }
        else {
          my %otherHit;
          $otherHit{'hitAccession'}    = $hit->accession;
          $otherHit{'hitSignificance'} = $hit->significance;
          $otherHit{'hitLength'}       = $hit->length;

          my $otherHits = $transcript{'otherHits'};
          push( @{$otherHits}, \%otherHit );
        }
      } # while hit
      push( @transcripts, \%transcript );
    }    # while results
    $filesDone++;
    print STDERR "\r".$filesDone."/".$fileCount;
  }    # for files
  print STDERR "\n";

  return \@transcripts;
}

__END__
