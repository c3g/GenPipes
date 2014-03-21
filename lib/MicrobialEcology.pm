#!/usr/bin/env perl

=head1 NAME

I<MicrobialEcology>

=head1 SYNOPSIS

MicrobialEcology->RDPWrapper()

=head1 DESCRIPTION

B<MicrobialEcology> This a library to analyze OTU data.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package MicrobialEcology;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Job;
use File::Basename;

# SUB
#-----------------------
sub RDPWrapper{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl'],
      ['rdp_classifier', 'moduleVersion.rdp_classifier']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' parallelRDP.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --RDP_training_set ' . LoadConfig::getParam($rH_cfg, 'DB', 'rdp_training_set', 1, 'filepath');
		$cmd .= ' --outfile ' . $outfile;
		$cmd .= ' --minWords ' . LoadConfig::getParam($rH_cfg, 'RDP', 'minWords', 1, 'int');
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'RDP', 'num_threads', 1, 'int');
		$cmd .= ' --fasta';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub RDPClassifier{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['rdp_classifier', 'moduleVersion.rdp_classifier'],
      ['java', 'moduleVersion.java']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' java -Xmx' . LoadConfig::getParam($rH_cfg, 'RDP', 'RAM', 1);
		$cmd .= ' rdp_classifier.jar';
		$cmd .= ' -q ' . $infile;
		$cmd .= ' -t ' . LoadConfig::getParam($rH_cfg, 'DB', 'rdp_training_set', 1, 'filepath');
		$cmd .= ' -o ' . $outfile;
		$cmd .= ' --minWords ' . LoadConfig::getParam($rH_cfg, 'RDP', 'minWords', 1, 'int');
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub addTaxToObs{
	my $rH_cfg             = shift;
	my $obs                = shift;
	my $rdp                = shift;
	my $outfile            = shift;
	my $outfileFailed      = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$obs] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
		$cmd .= 'module load '. LoadConfig::getParam($rH_cfg, 'memtime', 'moduleVersion.memtime').' && ';
		$cmd .= 'module load '.LoadConfig::getParam($rH_cfg, 'perl', 'moduleVersion.perl').' && ';
		$cmd .= 'module load ' . LoadConfig::getParam($rH_cfg, 'tools', 'moduleVersion.tools').' && ';	 
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['perl', 'moduleVersion.perl'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'addTaxToObs.pl';
		$cmd .= ' --seqobs ' . $obs;
		$cmd .= ' --rdp ' . $rdp;
		$cmd .= ' --cutoff ' . LoadConfig::getParam($rH_cfg, 'add_taxonomy', 'cutoff', 1, 'float');
		$cmd .= ' --outfile ' . $outfile;
		$cmd .= ' --outfile_failed ' . $outfileFailed;
		$cmd .= ' --tax_level ' . LoadConfig::getParam($rH_cfg, 'add_taxonomy', 'tax_level', 1);
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub filterOTUTable{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['perl', 'moduleVersion.perl'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' rmEmptyCol.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --outfile ' . $outfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub filterObsTable{
	my $rH_cfg             = shift;
	my $infileTsv          = shift;
	my $infileFasta        = shift;
	my $singletThreshold   = shift;
	my $outfileTsv         = shift;
	my $outfileFasta       = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infileTsv, $infileFasta] , [$outfileTsv, $outfileFasta]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['perl', 'moduleVersion.perl'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' rmEmptyCol.pl';
		$cmd .= ' --infile_tab ' . $infileTsv;
		$cmd .= ' --infile_fasta ' . $infileFasta;
		$cmd .= ' --outfile_tab ' . $outfileTsv;
		$cmd .= ' --outfile_fasta ' . $outfileFasta;
		$cmd .= ' --singlet_threhsold ' . $singletThreshold;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub splitOTUTable{
	my $rH_cfg            = shift;
	my $infile            = shift;
	my $matched           = shift;
	my $unmatched         = shift;
	my $selection         = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], [$matched]);

	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' splitOTUTable.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --matched ' . $matched;
		$cmd .= ' --unmatched ' . $unmatched;
		$cmd .= ' --select ' . $selection;	
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub convertOTUToBiom{
	my $rH_cfg             = shift;
	my $tsv                = shift;
	my $biom               = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$tsv] , [$biom]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime', 'moduleVersion.qiime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
    $cmd .= 'convert_biom.py ';
		$cmd .= ' -i ' . $tsv;
		$cmd .= ' -o ' . $biom;
		$cmd .= ' --biom_table_type=\\"otu table\\"';
		$cmd .= ' --process_obs_metadata taxonomy';
		$cmd .= ' && ';
		$cmd .= ' sed -e \'s/{/{\\n/g\' -e \'s/}/}\\n/g\' < '. $biom . ' > ' . $biom . '.tmp';
		$cmd .= ' && '; 
		$cmd .= ' mv ' . $biom . '.tmp ' . $biom;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub rarefy{
	my $rH_cfg             = shift;
	my $biom               = shift;
	my $biomRarefied       = shift;
	my $tsvRarefied        = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$biom] , [$biomRarefied, $tsvRarefied]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['perl', 'moduleVersion.perl'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'rarefaction.pl';
		$cmd .= ' --infile ' . $biom;
		$cmd .= ' --outfile ' . $biomRarefied; 
		$cmd .= ' --n ' . LoadConfig::getParam($rH_cfg, 'rarefaction', 'rarefactionThreshold', 1, 'int');
		$cmd .= ' && ';
		$cmd .=  'convert_biom.py ';
		$cmd .= ' -i ' . $biomRarefied; 
		$cmd .= ' -o ' . $tsvRarefied;
		$cmd .= ' -b --header_key=\\"taxonomy\\"';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub filterAndConvertToBiom{
	my $rH_cfg             = shift;
	my $tsv                = shift;
	my $tsvFiltered        = shift;
	my $frequency          = shift;
	my $threshold          = shift;
	my $biomFiltered       = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$tsv] , [$biomFiltered]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['perl', 'moduleVersion.perl'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'filterObsTable.pl';
		$cmd .= ' --otu_table ' . $tsv;
		$cmd .= ' --otu_table_out ' . $tsvFiltered;
		$cmd .= ' --frequency ' . $frequency;
		$cmd .= ' --threshold ' . $threshold;
		$cmd .= ' && ';
    $cmd .= 'convert_biom.py ';
		$cmd .= ' -i ' . $tsvFiltered;
		$cmd .= ' -o ' . $biomFiltered;
		$cmd .= ' --biom_table_type=\\"otu table\\"';
		$cmd .= ' --process_obs_metadata taxonomy';
		$cmd .= ' && ';
		$cmd .= ' sed -e \'s/{/{\\n/g\' -e \'s/}/}\\n/g\' < ' . $biomFiltered . ' > ' . $biomFiltered . '.tmp';
		$cmd .= ' && '; 
		$cmd .= ' mv ' . $biomFiltered . '.tmp ' . $biomFiltered;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub summarizeTaxonomyAbsolute{
	my $rH_cfg             = shift;
	my $biom               = shift;
	my $i                  = shift;
	my $outdir             = shift;
	
	my $basename = basename($biom);
	$basename =~ s{\.[^.]+$}{};

	my $ro_job = new Job();
	$ro_job->testInputOutputs([$biom], [$outdir."/".$basename."_L6.txt", $outdir."/".$basename."_L6.biom"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' summarize_taxa.py';
		$cmd .= ' -i ' . $biom;
		$cmd .= ' -L ' . $i;
		$cmd .= ' -o ' . $outdir;
		$cmd .= ' -a';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub summarizeTaxonomyRelative{
	my $rH_cfg             = shift;
	my $biom               = shift;
	my $i                  = shift;
	my $outdir             = shift;
	
	my $basename = basename($biom);
	$basename =~ s{\.[^.]+$}{};

	my $ro_job = new Job();
	$ro_job->testInputOutputs([$biom], [$outdir."/".$basename."_L6.txt", $outdir."/".$basename."_L6.biom"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' summarize_taxa.py';
		$cmd .= ' -i ' . $biom;
		$cmd .= ' -L ' . $i;
		$cmd .= ' -o ' . $outdir;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub plotTaxa{
	my $rH_cfg             = shift;
	my $infiles            = shift;
	my $outdir             = shift;
	
	my @infiles = split(/,/, $infiles);

	my $ro_job = new Job();
	$ro_job->testInputOutputs(\@infiles, ["$outdir/bar_charts.html"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'plot_taxa_summary.py';
		$cmd .= ' -i ' . $infiles;
		$cmd .= ' -o ' . $outdir;
		$cmd .= ' -c bar';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub calculateAbundanceThresholds{
	my $rH_cfg             = shift;
	my $tsv                = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$tsv] , [""]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['perl', 'moduleVersion.perl'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'abundanceThreshold.pl';
		$cmd .= ' --infile ' . $tsv;
		$cmd .= ' --outdir ' . $outdir;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub phylumBarplot{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outdir             = shift;
	my $prefix             = shift;
	#my $outfileGraph       = shift;
	#my $outfileTable       = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], [$outdir."/".$prefix.".pdf"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['R', 'moduleVersion.R'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' phylumBarplot.R';
		$cmd .= ' -i ' . $infile;
		$cmd .= ' -o ' . $outdir;
		$cmd .= ' -p ' . $prefix;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub filterForPyNAST{
	my $rH_cfg             = shift;
	my $tsv                = shift;
	my $fasta              = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$tsv, $fasta] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['perl', 'moduleVersion.perl'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' filterForPyNAST.pl';
		$cmd .= ' --infile_otu_table ' . $tsv;
		$cmd .= ' --infile_fasta ' . $fasta;
		$cmd .= ' > ' . $outfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub PyNAST{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $log                = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['python', 'moduleVersion.python'],
      ['pynast', 'moduleVersion.pynast'],
      ['openmpi', 'moduleVersion.openmpi']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' mpirun -np ' . LoadConfig::getParam($rH_cfg, 'pynast', 'num_threads', 1, 'int');
		$cmd .= ' pynast ';
		$cmd .= ' -i ' . $infile; 
		$cmd .= ' -p 10';
		$cmd .= ' -l 50';
		$cmd .= ' -g ' . $log;
		$cmd .= ' -a ' . $outfile;
		$cmd .= ' -t ' . LoadConfig::getParam($rH_cfg, 'DB',  'core', 1, 'filepath'); 
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub filterPyNASTAlignment{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outdir."/pynast_pfiltered.fasta"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' filter_alignment.py';
		$cmd .= ' -i ' . $infile;
		$cmd .= ' -o ' . $outdir;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub fasttree{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['fasttree', 'moduleVersion.fasttree']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'FastTree';
		$cmd .= ' -nt ' . $infile . ' > '. $outfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub betaDiversity{
	my $rH_cfg             = shift;
	my $biom               = shift;
	my $distance           = shift;
	my $outdir             = shift;
	my $tree               = shift; 
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$biom], [$outdir."/"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
    $cmd .= 'unset LD_LIBRARY_PATH; ';
		$cmd .=	' memtime ';
		$cmd .= 'beta_diversity.py';
		$cmd .= ' -i ' . $biom;
		$cmd .= ' -m ' . $distance;
		$cmd .= ' -o ' . $outdir;
		$cmd .= ' -t ' . $tree;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub principalCoordinates{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'principal_coordinates.py';
		$cmd .= ' -i ' . $infile;
		$cmd .= ' -o ' . $outfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub pca2dPlot{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outdir             = shift;
	
	my $basename = basename($infile);
	$basename =~ s{\.[^.]+$}{};
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outdir."/".$basename."_2D_PCoA_plots.html"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'make_2d_plots.py';		
		$cmd .= ' -i ' . $infile;
		$cmd .= ' -m ' . LoadConfig::getParam($rH_cfg, 'default', 'mappingFile', 1, 'filepath');
		$cmd .= ' -o ' . $outdir;

		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub pca3dPlot{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outdir             = shift;
	
	my $basename = basename($infile);
	$basename =~ s{\.[^.]+$}{};
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outdir."/".$basename."_3D_PCoA_plots.html"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'make_3d_plots.py';		
		$cmd .= ' -i ' . $infile;
		$cmd .= ' -m ' . LoadConfig::getParam($rH_cfg, 'default', 'mappingFile', 1, 'filepath');
		$cmd .= ' -o ' . $outdir;

		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub multipleRarefaction{
	my $rH_cfg             = shift;
	my $biom               = shift;
	my $tsv                = shift;
	my $outdir             = shift;
	my $rootOutdir         = shift;
	
	my $ro_job = new Job();
  my $dummyOutfile = "$rootOutdir/rarefactions.done";
	$ro_job->testInputOutputs([$biom] , [$dummyOutfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['perl', 'moduleVersion.perl'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .= ' rm -rf ' . $rootOutdir . '/rarefaction/* && ';
		$cmd .= ' rm -rf ' . $rootOutdir . '/alpha_rarefaction/* && ';
		$cmd .= ' rm -rf ' . $rootOutdir . '/collated/* && ';
		$cmd .= ' rm -rf ' . $rootOutdir . '/plots/* && ';
		$cmd .=	' memtime ';
		$cmd .= 'multipleRarefactions.pl';
		$cmd .= ' --infile ' . $biom;
		$cmd .= ' --infileTsv ' . $tsv;
		$cmd .= ' --outdir ' . $outdir;
		$cmd .= ' --threshold ' . LoadConfig::getParam($rH_cfg, 'multiple_rarefaction', 'minFractionThreshold', 1, 'float');
		$cmd .= ' --n ' . LoadConfig::getParam($rH_cfg, 'multiple_rarefaction', 'rarefactionThreshold', 1, 'int');
		$cmd .= ' --step ' . LoadConfig::getParam($rH_cfg, 'multiple_rarefaction', 'step', 1, 'int');
		$cmd .= ' --permutations ' . LoadConfig::getParam($rH_cfg, 'multiple_rarefaction', 'perm', 1, 'int');
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'multiple_rarefaction', 'num_threads', 1, 'int');
    $cmd .= ' && touch ' . $dummyOutfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub alphaDiversity{
	my $rH_cfg             = shift;
	my $indir              = shift;
	my $outdir             = shift;
  
  my $dummyOutfile = "$outdir/../alphaDiv.done";
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([""] , [$dummyOutfile]);
	
	my $m = join(',', LoadConfig::getParam($rH_cfg, 'alpha_diversity', 'm'));
	$m =~ s/:/,/g;

	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
    $cmd .= 'unset LD_LIBRARY_PATH; ';
		$cmd .=	' memtime ';
		$cmd .= 'parallel_alpha_diversity.py';
		$cmd .= ' -i ' . $indir;
		$cmd .= ' -o ' . $outdir;
		$cmd .= ' -m ' . $m;
		$cmd .= ' -O ' . LoadConfig::getParam($rH_cfg, 'alpha_diversity', 'num_threads', 1, 'int');
    $cmd .= ' && touch ' . $dummyOutfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub collateAlpha{
	my $rH_cfg             = shift;
	my $indir              = shift;
	my $outdir             = shift;
  
  my $dummyOutfile = "$outdir/../collateAlpha.done";
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([""] , [$dummyOutfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'collate_alpha.py';
		$cmd .= " -i ".$indir;
		$cmd .= " -o ".$outdir;
    $cmd .= ' && touch ' . $dummyOutfile;
		
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub rarefactionPlots{
	my $rH_cfg             = shift;
	my $indir              = shift;
	my $outdir             = shift;
	
  my $dummyOutfile = "$outdir/rarefactionPlots.done";
	
  my $ro_job = new Job();
	$ro_job->testInputOutputs([""] , [$dummyOutfile, "$outdir/rarefactionPlots.pdf"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python'],
      ['R', 'moduleVersion.R'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'make_rarefaction_plots.py';
		$cmd .= ' -i ' . $indir;
		$cmd .= ' -m ' . LoadConfig::getParam($rH_cfg, 'default', 'mappingFile', 1, 'filepath');
		$cmd .= ' -o ' . $outdir;
    $cmd .= ' && ';
    $cmd .= 'rarefactionPlots.R';
    $cmd .= ' -r ' . $outdir . '/average_tables/observed_speciesSampleID.txt';
    $cmd .= ' -o ' . $outdir;
    $cmd .= ' -p rarefactionPlots';
    $cmd .= ' && touch ' . $dummyOutfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub upgmaClustering{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfile            = shift;
	
  my $outdir = dirname $outfile;
  my $prefix = basename $outfile;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['qiime-dependencies', 'moduleVersion.qiime-dependencies'],
      ['qiime', 'moduleVersion.qiime'],
      ['python', 'moduleVersion.python'],
      ['R', 'moduleVersion.R'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'upgma_cluster.py';
		$cmd .= ' -i ' . $infile;
		$cmd .= ' -o ' . $outfile;
    $cmd .= ' && memtime ';
    $cmd .= 'plotTree.R';
    $cmd .= ' -t ' . $outfile;
    $cmd .= ' -o ' . $outdir;
    $cmd .= ' -p ' . $prefix;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub otuHeatMap{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outdir             = shift;
  my $prefix             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , ["$outdir/$prefix.pdf"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['R', 'moduleVersion.R'],
      ['tools', 'moduleVersion.tools']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'OTUheatmap.R';
		$cmd .= ' -t ' . $infile;
    $cmd .= ' -o ' . $outdir;
		$cmd .= ' -p ' . $prefix;
		$cmd .= ' -n ' . LoadConfig::getParam($rH_cfg, 'otu_heatmap', 'n', 1, 'int');
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub template4{
	my $rH_cfg             = shift;
	my $tsv                = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs(undef , undef);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime']
    ]) . ' && ';
		$cmd .=	' memtime ';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub template1{
	my $rH_cfg             = shift;
	my $tsv                = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs(undef , undef);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime']
    ]) . ' && ';
		$cmd .=	' memtime ';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}
1
