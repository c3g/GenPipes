#!/usr/env/perl

=head1 NAME

I<HOMER>

=head1 SYNOPSIS

HOMER->parseGenome()
HOMER->makeTagDirectory()
HOMER->findMotifs()
HOMER->annotatePeaks()
HOMER->makeUCSCFile ()

=head1 DESCRIPTION

B<HOMER> is a software for motif discovery and next-gen sequecing analysis

Input = file_name

Output = array

=head1 AUTHOR

B<Johanna Sandoval> - I<johanna.sandoval@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package HOMER;

# Strict Pragmas
#--------------------------
use strict;
use warnings;
use LoadConfig;
#use HomerConfig;
#--------------------------

# Dependencies
#-----------------------

# SUB
#-----------------------

sub parseGenome {
  my $rH_cfg     = shift ;
  my $genomeName ;
  # Retrieve genome name from ini file
  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'genomeName');
  
  # Check if genome exists in homer config file. 
	my $config = HomerConfig::loadConfigFile();
	
  if (!exists({$config->{GENOMES}}->{$refGenome} )) {
    print STDERR "\n#WARNING: Genome $refGenome not found in Homer config.txt file \n#QC, annotations and Motif analysis will not be executed\n\n";
    $genomeName = 'none';
	} else {
		$genomeName = $refGenome;
	}
	return $genomeName;
}

sub makeTagDirectory {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $sortedBAM       = shift;
  my $outputDir       = shift;
  my $processUnmapped = shift;

  my $refGenome = LoadConfig::getParam($rH_cfg, 'default', 'genomeName');
  #my $refGenome = parseGenome( $rH_cfg );
  my @commands;
  my $command;
  if ( defined $refGenome ) {
		
		my $tmpSampleFile = LoadConfig::getParam($rH_cfg, 'qcTags', 'tmpDir').'/dmp-XXXXXX.'.$sampleName;
		my $tmpSampleFileName = LoadConfig::getParam($rH_cfg, 'qcTags', 'tmpDir').'/filename.'.$sampleName;
 		#$command .= ' rm  -f '.$tmpSampleFile.' ;'; 
 		$command .= ' tmpFile=\$(mktemp -u '.$tmpSampleFile.'); echo \$tmpFile > '.$tmpSampleFileName.';'; 
		push (@commands, $command);
 		$command = ' module load '.LoadConfig::getParam($rH_cfg, 'qcTags', 'moduleVersion.samtools').' ;'; 
 		$command .= ' tmpFile=\$(cat '.$tmpSampleFileName.'); '; 
 		$command .= ' mkfifo \$tmpFile && '.' samtools view -h '.$sortedBAM.' > \$tmpFile ;'; 
 		push (@commands, $command);
 		$command  = ' module load '. LoadConfig::getParam($rH_cfg, 'qcTags', 'moduleVersion.python').' '.LoadConfig::getParam($rH_cfg, 'qcTags', 'moduleVersion.homer'). ';';
 		$command .= ' tmpFile=\$(cat '.$tmpSampleFileName.'); '; 
 		$command .= ' makeTagDirectory '.$outputDir.'/'. $sampleName.' \$tmpFile -checkGC -genome '.$refGenome.' ; ';
 		$command .= ' rm  -f '.$tmpSampleFile.'; ';
 		push (@commands, $command);
# 		$command .= ' module load '.LoadConfig::getParam($rH_cfg, 'qcTags', 'moduleVersion.samtools').' '. LoadConfig::getParam($rH_cfg, 'qcTags', 'moduleVersion.python').' '.LoadConfig::getParam($rH_cfg, 'qcTags', 'moduleVersion.homer'). ';'; 
# 		$command .= ' rm  -f '.$tmpSampleFile.' ;'; 
# 		$command .= ' mktemp -u '.$tmpSampleFile.' && mknod '.$tmpSampleFile.' p && '.' samtools view -h '.$sortedBAM.' > '.$tmpSampleFile.' & '; 
# 		$command .= ' makeTagDirectory '.$outputDir.'/'. $sampleName.' '.$tmpSampleFile.' -checkGC -genome '.$refGenome.' ;';		
# 		push (@commands, $command);

  }else{
		@commands=();
		print STDERR "\n#WARNING: Genome $refGenome not defined \n#QC, annotations and Motif analysis will not be executed\n\n";
  }  
  return \@commands;
}


sub makeUCSCFile {
	my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $tagDirectory    = shift;
  my $outputWiggle    = shift;
  
	my $command = undef;
	
  my $latestOutputFile = -M $outputWiggle;
  
  if(!defined($latestOutputFile ))
  {
	  $command .= ' module load ' .LoadConfig::getParam($rH_cfg, 'default' , 'moduleVersion.python') .' ' . LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.homer'). ';';
	  $command .= ' makeUCSCfile ' . $tagDirectory .' -o ' . $outputWiggle .';';
	}
	
	return $command;  
}

sub annotatePeaks{
	my $rH_cfg        = shift;
	my $designName    = shift;
	my $InputBed   		 = shift;
  my $outputDir     = shift;
  my $command ;
  my $genomeName    = LoadConfig::getParam($rH_cfg, 'annotation', 'genomeName');
  
  
	$command .= ' module load ' .LoadConfig::getParam($rH_cfg, 'default' , 'moduleVersion.python') .' ' . LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.homer'). ';';
	$command .= ' annotatePeaks.pl '.$InputBed.' '. $genomeName . ' -gsize '.$genomeName.' -cons -CpG -go '.$outputDir.'/'.$designName.' -genomeOntology ' .$outputDir.'/'.$designName. ' > '.$outputDir.'/'.$designName.'.annotated.csv'.';';
	return $command;  
  
}

 sub generateMotif{
	my $rH_cfg        = shift;
	my $designName    = shift;
	my $InputBed   		 = shift;
  my $outputDir     = shift;
  my $command ;
  my $genomeName    = LoadConfig::getParam($rH_cfg, 'motif', 'genomeName');
  
  my $optionsThreads ;
  if (defined(LoadConfig::getParam( $rH_cfg, 'motif', 'homermotifThreads')) &&  LoadConfig::getParam( $rH_cfg, 'motif', 'homermotifThreads') ne "" ){
			$optionsThreads = '-p '.LoadConfig::getParam( $rH_cfg, 'motif', 'homermotifThreads').' ';
  }else{
			$optionsThreads = ' ';
  }
	$command .= ' module load ' .LoadConfig::getParam($rH_cfg, 'default' , 'moduleVersion.python') .' ' . LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.homer'). ' '. LoadConfig::getParam($rH_cfg, 'default', 'moduleVersion.weblogo').';';
	$command .= ' findMotifsGenome.pl '.$InputBed.' '. $genomeName . ' ' .$outputDir. ' ' . $optionsThreads.';';
	return $command;  
  
}

1;
