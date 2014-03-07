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

    # Check modules availability
    checkModules(\%config);

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

# Check by a system call if all modules defined in config file are available
sub checkModules {
  my $rH_cfg = shift;

  # Retrieve all module values in config file assuming that all module key names contain "moduleVersion"
  for my $key (keys %$rH_cfg) {
    if ($key =~ /moduleVersion/) {
      my $moduleValue = $rH_cfg->{$key};

      # Bash "module" cmd not found when switching to a different shell so source module.sh is required
      my $moduleShowOutput = `bash -c 'source /etc/profile.d/modules.sh; module show $moduleValue 2>&1'`;
      $moduleShowOutput !~ /Error/i or die "Error in config file with $moduleValue:\n$moduleShowOutput";
    }
  }
}

# Retrive param in config file with optional definition check and type validation
sub getParam {
  my $rH_cfg  = shift;
  my $section = shift;
  my $value = shift;
  my $required = shift;
  my $type = shift;

  unless (defined($required)) {$required = 1}; # By default, parameter is required to be defined in the config file

  # Search param in specified section
  my $retVal = $rH_cfg->{$section . '.' . $value};

  # If not found, search param in default section
  unless (defined($retVal)) {$retVal = $rH_cfg->{'default.' . $value}};

  if (defined($retVal)) {
    if (defined($type)) {
      if ($type eq 'int') {
        $retVal =~ /^[-+]?\d+$/ or die "Error: parameter \"[" . $section . "] " . $value . "\" value $retVal is not an integer!";
      } elsif ($type eq 'float') {
        $retVal =~ /^[-+]?\d*\.?\d+([eE][-+]?\d+)?$/ or die "Error: parameter \"[" . $section . "] " . $value . "\" value $retVal is not a float!";
      } elsif ($type eq 'path') {
        my $tmpVal = `echo -n $retVal`;
        -e $tmpVal or die "Error: parameter \"[" . $section . "] " . $value . "\" value $retVal is not a valid path!";
      } elsif ($type eq 'filepath') {
        my $tmpVal = `echo -n $retVal`;
        -f $tmpVal or die "Error: parameter \"[" . $section . "] " . $value . "\" value $retVal is not a valid plain file path!";
      } elsif ($type eq 'dirpath') {
        my $tmpVal = `echo -n $retVal`;
        -d $tmpVal or die "Error: parameter \"[" . $section . "] " . $value . "\" value $retVal is not a valid directory path!";
      } elsif ($type eq 'array') {
        if (ref($retVal) ne "ARRAY") {
          if (ref(\$retVal) eq "SCALAR") {
            # Convert scalar into 1-element array due to "tie" function ambiguity for 1-value list without ","
            $retVal = [$retVal];
          } else {
            die "Error: parameter \"[" . $section . "] " . $value . "\" value $retVal is not an array!";
          }
        }
      } else {
        die "Error: parameter \"[" . $section . "] " . $value . "\" validation type \"$type\" unknown!";
      }
    }
  } elsif ($required) {
    die "Error: parameter \"[" . $section . "] " . $value . "\" is not defined in config file!";
  } else {
    $retVal = "";
  }

  return $retVal;
}

# Return "module load ..." command string from the given list of modules as [[section, moduleVersion], ...]
sub moduleLoad {
  my $rH_cfg = shift;
  my $rA_modules = shift;

  # Retrieve module values from the configuration file
  my @moduleValues = map {getParam($rH_cfg, $_->[0], $_->[1], 1)} @$rA_modules;

  return "module load " . join(" ", @moduleValues);
}

1;
