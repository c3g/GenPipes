#!/usr/env/perl

=head1 NAME

I<SequenceDictionaryParser>

=head1 SYNOPSIS

SequenceDictionaryParser->readDictFile(Config_file_name)

=head1 DESCRIPTION

B<SequenceDictionaryParser> is a library that reads from a sequence dictionary file and 
returns the globlal values in a hash.

=head1 AUTHOR

B<Louis Letourneau> - I<louis.letourneau@mail.mcgill.ca>

=head1 DEPENDENCY

=cut

package SequenceDictionaryParser;

# Strict Pragmas
#---------------------
use strict;
use warnings;

#---------------------

# Dependencies
#--------------------
use LoadConfig;

#--------------------

# SUB
#--------------------
sub readDictFile {
  my $rH_cfg = shift;
  my %dictionary;

  my $refDictFile = LoadConfig::getParam( $rH_cfg, 'default', 'referenceSequenceDictionary' );

  if ( -e $refDictFile ) {
    open(FILE, $refDictFile) or die "Cannot open ".$refDictFile."\n";
    while(my $line = <FILE>) {
      if($line =~ /^\@SQ\tSN:([^\t]+)\tLN:(\d\+)/) {
        $dictionary{$1} = {'name' => $1, 'size' => $2};
      }
    }
    close(FILE);
    return \%dictionary;

  }
  else {
    print "Reference dictionary absent. First create one with picard\n";
    exit 1;
  }
}

1;
