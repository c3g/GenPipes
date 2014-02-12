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
  my ($self, $file) = @_;
  my %config;

  # Load in the local config for the invoking script if it exists.
  #---------------------------------------------------------------
  if (-e "$file") {
    tie %config, "Config::Simple", $file;
    tied(%config)->autosave(1);
    return %config;
  }

  # Prompt an error message othewise
  #---------------------------------
  else {
    die "Error: config file $file does not exit or is not a valid file!\n\n" .
      "First, create a config file according to this model:\n\n" .
      "[section]\nNAME= Value\n\nEXAMPLE\n\n[PATH]\nJAVA=/opt/jdk/jdk1.6.0_33/bin/java\n\n" .
      "In this case a hash \n\$hash{PATH.JAVA}=/opt/jdk/jdk1.6.0_33/bin/java\n\nwill be returned\n";
  }
}

sub getParam {
  my $rH_cfg  = shift;
  my $section = shift;
  my $value   = shift;
  my $validation = shift; # If set to true, raise an error if parameter is not defined in the config file (default: false)

  my $retVal = $rH_cfg->{$section . '.' . $value};
  if (!defined($retVal)) {
    $retVal = $rH_cfg->{'default.' . $value};
    if (!defined($retVal)) {
      if ($validation) {
        die "Error: parameter \"[" . $section . "] " . $value . "\" is not defined in config file!";
      } else {
        $retVal = "";
      }
    }

    #    if(ref($retVal) eq "ARRAY" && scalar(@{$retVal}) == 0) {
    #      $retVal = undef;
    #    }
  }
  return $retVal;
}

# Return "module load ..." command string from the given list of modules as [[section, moduleVersion], ...]
sub moduleLoad {
  my $rH_cfg = shift;
  my $rA_modules = shift;

  # Retrieve module values from the configuration file
  my @moduleValues = map {getParam($rH_cfg, $_->[0], $_->[1], 1)} @$rA_modules;

  # Check by a system call if module is available
  for my $moduleValue (@moduleValues) {
    my $moduleShowOutput = `source /etc/profile.d/modules.sh; module show $moduleValue 2>&1`;
    $moduleShowOutput !~ /Error/i or die "Error in config file:\n$moduleShowOutput";
  }

  return "module load " . join(" ", @moduleValues);
}

1;
