#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Iterator::FastaDb;
use Iterator::FastqDb;
use Data::Dumper;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
hpss_itags.pl

PURPOSE:
To loop recursively through a directory and compress each .fastq
to .gz file. 

INPUT:
--indir <string>  : Directory where itags runs are to be backed-up.

OUTPUT:
No output options.

NOTES:

BUGS/LIMITATIONS:
 
AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

ENDHERE

## OPTIONS
my ($help, $indir);
my $verbose = 0;

GetOptions(
    'indir=s' 	=> \$indir,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

## Validate
die "--indir arg missing\n" unless($indir);

## It is intended to recursively through a directory tree
## and compress all *.fastq in *.fastq.gz


## MAIN
sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 
		
		if(substr($filename, -6) eq ".fastq"){
			print STDOUT "Compressing ".$fullpath." into .gz archive...\n";
			system("gzip ".$fullpath);
			$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
		}
	}
}

## MAIN

# Compress .fastq into .gz
find (\&eachFile, $indir);

exit;
