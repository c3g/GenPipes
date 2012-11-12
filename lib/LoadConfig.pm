#!/usr/env/perl

=head1 NAME

I<LoadConfig>

=head1 SYNOPSIS

LoadConfig->readConfigFile(Config_file_name)

=head1 DESCRIPTION

B<LoadConfig> is a library that reads from a config file and 
returns the globlal values in a hash.

Input = /path/Config_file_name

Output = %hash 


The configuration file must be divided into section.
Each of which will hold a key = value option.

B<For example:>


#########################################

[path]

BlastDb = /path/to/Blastdb

[file]

Fastq = /path/to/fataq/file.gz

#####################################


These fields would be exported as:

$config{'path.BlastDb'} = /path/to/Blastdb

$config{'file.Fastq'} = /path/to/fataq/file.gz


=head1 AUTHOR

B<David Morais> - I<dmorais@cs.bris.ac.uk>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

B<Config::Simple> Used to parse config file

B<File::Basename> path parsing

B<Cwd> path parsing

=cut

package LoadConfig;

# Strict Pragmas
#---------------------
use strict;
use warnings;
#---------------------

# Dependencies
#--------------------
use Data::Dumper;
use Config::Simple;
use File::Basename;
use Cwd 'abs_path';
#--------------------


# SUB
#--------------------
sub readConfigFile {
    my ( $self, $file ) = @_;
    my %config;

	
    #Load in the local config for the invoking script if it exists.
    #---------------------------------------------------------------
    if ( -e "$file" ) {
        tie %config, "Config::Simple", $file;
        tied(%config)->autosave(1);
        return %config;

    }

	
    #Prompt an error message othewise
    #---------------------------------
    else {
        print "Firts, create a config file in the same \n" .
          "The file must be created according to this model\n\n" .
          "[section]\nNAME= Value\n\nEXAMPLE\n\n[PATH]\nJAVA=/opt/jdk/jdk1.6.0_33/bin/java\n\n" .
          "In this case a hash  \n\$hash{PATH.JAVA}=/opt/jdk/jdk1.6.0_33/bin/java\n\nwill be returned\n";
        exit;

    }

}

1;
