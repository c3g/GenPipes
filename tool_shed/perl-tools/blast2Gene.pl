#!/usr/bin/perl -w

use strict;
use Bio::SearchIO;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::DB::EUtilities;

&main();

sub main {
  my %orfs;
  parseGlimmer( \%orfs );
  blast( \%orfs );
  getGenes( \%orfs );

  for my $rH_orf ( values %orfs ) {
    my $glimmerName      = $rH_orf->{'name'};
    my $start            = $rH_orf->{'start'};
    my $end              = $rH_orf->{'stop'};
    my $score            = $rH_orf->{'glimmerScore'};
    my $blastEValue      = $rH_orf->{'hitSignificance'};
    my $blastScore       = $rH_orf->{'hitScore'};
    my $geneDescription  = $rH_orf->{'hitDescription'};
    my $geneName         = $rH_orf->{'gene'}->{'name'};
    my $geneID           = $rH_orf->{'gene'}->{'id'};
    my $proteinAccession = $rH_orf->{'hitAccession'};

    if ( $end < $start ) {
      print("     CDS              complement($start..$end)\n");
    }
    else {
      print("     CDS              $start..$end\n");
    }
    print("                      \/note=\"predicted using Glimmer.\n");
    if ( defined $blastEValue ) {
      print("                      $geneDescription\n");
      print("                      BlastEValue: $blastEValue\n");
      print("                      BlastScore: $blastScore\n");
    }
    print("                      Glimmer Score: $score \n");
    print("                      \"\n");
    print("                      \/transl_table=11\n");
    print("                      \/codon_start=1\n");
    print("                      \/locus_tag=\"$glimmerName\"\n");
    if ( defined $proteinAccession ) {
      print("                      \/protein_id=\"$proteinAccession\"\n");
      if ( defined $geneName ) {
        print("                      \/gene=\"$geneName\"\n");
        print("                      \/db_xref=\"GeneID:$geneID\"\n");
      }
    }
  }
}

sub parseGlimmer() {
  my $rHoH_orf = shift;

  my $glimmerOutput = $ARGV[0];
  open( GLIM, "<$glimmerOutput" ) or die "Can't open glimmer file\n";
  while (<GLIM>) {
    chomp;
    if ( substr( $_, 0, 1 ) eq '>' ) {
      next;
    }

    my @results = split;
    my %orf;
    $orf{'name'}               = $results[0];
    $orf{'start'}              = $results[1];
    $orf{'stop'}               = $results[2];
    $orf{'glimmerScore'}       = $results[4];
    $rHoH_orf->{ $results[0] } = \%orf;
  }
  close(GLIM);
}

sub blast {
  my $rHoH_orf = shift;

  for ( my $blastFileIdx = 1 ; $blastFileIdx < @ARGV ; $blastFileIdx++ ) {
    my $blastOutput = $ARGV[$blastFileIdx] or die "Can't open blast file\n";

    # comment out the next line to read STDIN
    my $in = new Bio::SearchIO(
      -format => 'blastxml',
      -file   => $blastOutput
    );

    while ( my $result = $in->next_result ) {
      my @query = split( / /, $result->query_description, );
      my $queryName = substr( $query[0], 0, length( $query[0] ) - 2 );

      my $rH_orf = $rHoH_orf->{$queryName};

      my @stats  = $result->available_statistics;
      my @params = $result->available_parameters;

      # Use only the first, best, hit.
      my $hit = $result->next_hit;
      if ( defined $hit ) {
        $hit->name =~ /gi\|(\d+)/;
        $rH_orf->{'hitGI'}              = $1;
        $rH_orf->{'hitAccession'}       = $hit->accession;
        $rH_orf->{'hitFullDescription'} = $hit->description;
        $rH_orf->{'hitSignificance'}    = $hit->significance;
        $rH_orf->{'hitScore'}           = $hit->raw_score;

        my @entries = split( /\|/, $rH_orf->{'hitFullDescription'} );
        $rH_orf->{'hitDescription'} = $entries[0];
      }
    }
  }
}

sub getGenes {
  my $rHoH_orf = shift;

  my %seenGIs;
  for my $rH_orf ( values %$rHoH_orf ) {
    if ( !defined $rH_orf->{'hitGI'} ) {
      next;
    }
    $seenGIs{ $rH_orf->{'hitGI'} } = 1;
  }
  my @allGIs = keys %seenGIs;

  my $startIdx = 0;
  my $endIdx   = 0;
  while (1) {
    if ( $startIdx > @allGIs ) {
      last;
    }
    $endIdx += 100;
    if ( $endIdx > @allGIs ) {
      $endIdx = @allGIs;
    }
    my @gis = @allGIs[ $startIdx .. $endIdx ];
    $startIdx += 100;
    my $eutil = Bio::DB::EUtilities->new(
      -eutil          => 'elink',
      -email          => 'louis.letourneau@mail.mcgill.ca',
      -db             => 'gene',
      -verbose        => 1,
      -dbfrom         => 'protein',
      -id             => \@gis,
      -correspondence => 1
    );

    my %proteins;
    my %genes;

    # iterate through the LinkSet objects
    while ( my $ds = $eutil->next_LinkSet ) {
      my @retIDs = $ds->get_ids();
      if ( @retIDs > 0 ) {
        if ( @retIDs > 1 ) {
          print "NbGenes for ("
            . join( ',', $ds->get_submitted_ids() ) . "): "
            . @retIDs . "\n";
        }
        my %gene;
        $gene{'id'} = $retIDs[0];

        $genes{ $gene{'id'} } = \%gene;
        for my $subIDs ( $ds->get_submitted_ids() ) {
          $proteins{$subIDs} = \%gene;
        }
      }
      else {
        print "No genes found for: "
          . join( ',', $ds->get_submitted_ids ) . "\n";
      }
    }

    $eutil->reset_parameters(
      -eutil   => 'elink',
      -email   => 'louis.letourneau@mail.mcgill.ca',
      -db      => 'gene',
      -verbose => 1,
      -cmd     => 'neighbor_history',
      -dbfrom  => 'protein',
      -id      => \@gis,
    );
    my $hist = $eutil->next_History || die "No history data returned";

    $eutil->set_parameters(
      -eutil   => 'esummary',
      -history => $hist
    );

    # iterate through the individual DocSum objects (one per ID)
    while ( my $ds = $eutil->next_DocSum ) {
      my $rH_gene = $genes{ $ds->get_id() };

      # flattened mode, iterates through all Item objects
      while ( my $item = $ds->next_Item('flattened') ) {
        if ( $item->get_content ) {
          if ( $item->get_name eq 'Name' ) {
            $rH_gene->{'name'} = $item->get_content();
          }
          elsif ( $item->get_name eq 'Description' ) {
            $rH_gene->{'description'} = $item->get_content();
          }
        }
      }
    }

    for my $rH_orf ( values %$rHoH_orf ) {
      if ( !defined $rH_orf->{'hitGI'} ) {
        next;
      }
      $rH_orf->{'gene'} = $proteins{ $rH_orf->{'hitGI'} };
    }
  }
}

__END__
