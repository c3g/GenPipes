
#!/usr/bin/env perl

=head1 NAME

I<RRNAAMplicons>

=head1 SYNOPSIS

RRNAAmplicons->dnaclust()

=head1 DESCRIPTION

B<RRNAAmplicons> This a library to analyze rRNA amplicons using the RRNAAMplicons pipeline.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package RRNAAmplicons;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use Job;

# SUB
#-----------------------
sub mergeBarcodes{
	my $rH_cfg             = shift;
	my $rHoAoH_sampleInfo  = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([""], [$outdir.'/raw_reads/paired.fastq']);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
		my $cmd_R1 = '';
		my $cmd_R2 = '';

		# Fofn
		# Concatenate all libraries. Loop through .fastq.gz that will be specified by sampleSetup.pl.
		# Basically, concatenate all *R1.fastq.gz and *R2.fastq.gz separately. Assuming barcodes
		# are in headers.
		# Loop through each samples in sample sheet
		#$cmd .= 'module load '. LoadConfig::getParam($rH_cfg, 'memtime', 'moduleVersion.memtime').' && ';
		my $rawReadDir  = LoadConfig::getParam($rH_cfg, 'default','rawReadDir', 1, 'dirpath');
    $cmd_R1 .= LoadConfig::moduleLoad($rH_cfg, [
      ['perl', 'moduleVersion.perl'],
      ['tools', 'moduleVersion.tools'],
      ['memtime', 'moduleVersion.memtime']
    ]) . ' && memtime cat ';
		$cmd_R2 .= " memtime cat ";
		for my $sampleName (keys %{$rHoAoH_sampleInfo}) {
			my $rAoH_sampleLanes = $rHoAoH_sampleInfo->{$sampleName};
			for my $rH_laneInfo (@$rAoH_sampleLanes) { #this rH_laneInfo contains the complete line info from the sample sheet for this sample.
				$cmd_R1 .= " " . $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' . $sampleName . "." . $rH_laneInfo->{'libraryBarcode'} . "." . LoadConfig::getParam($rH_cfg, 'default','qual', 1, 'int') . ".pair1.fastq.gz";
				$cmd_R2 .= " " . $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' . $sampleName . "." . $rH_laneInfo->{'libraryBarcode'} . "." . LoadConfig::getParam($rH_cfg, 'default','qual', 1, 'int') . ".pair2.fastq.gz";
				#$cmd_R1 .= " " . $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' . $sampleName . "." . LoadConfig::getParam($rH_cfg, 'default','qual', 1, 'int') . ".pair1.fastq.gz";
				#$cmd_R2 .= " " . $rawReadDir .'/' .$sampleName .'/run' .$rH_laneInfo->{'runId'} . "_" . $rH_laneInfo->{'lane'} .'/' . $sampleName . "." . LoadConfig::getParam($rH_cfg, 'default','qual', 1, 'int') . ".pair2.fastq.gz";
			}
		}
		$cmd_R1 .= " | gunzip -c > ".$outdir."/raw_reads/reads_1.fastq";
		$cmd_R2 .= " | gunzip -c > ".$outdir."/raw_reads/reads_2.fastq";
		
		$cmd = $cmd_R1 . ' && ' . $cmd_R2 . ' && ';

		# Then merge as an interleaved file. 
		$cmd .= " memtime ";
		$cmd .= " joinPairedFastqs.pl";
		$cmd .= " --reads_1 ".$outdir."/raw_reads/reads_1.fastq";
		$cmd .= " --reads_2 ".$outdir."/raw_reads/reads_2.fastq";
		$cmd .= " --outfile ".$outdir."/raw_reads/paired.fastq";
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub dukWrapper{
	my $rH_cfg           = shift;
	my $infileFastq      = shift;
	my $contam           = shift;
	my $ncontam          = shift;
	my $log              = shift;
	my $db               = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infileFastq] , [$contam, $ncontam]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl'],
      ['duk', 'moduleVersion.duk']
    ]) . ' && ';
		$cmd .= ' memtime ';	 
		$cmd .= ' contamWrapper.pl';
		$cmd .= ' --infile ' . $infileFastq;
		$cmd .= ' --outfile_matched ' . $contam;
		$cmd .= ' --outfile_unmatched ' . $ncontam;
		$cmd .= ' --log ' . $log;
		$cmd .= ' --db ' . $db;
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'duk_wrapper', 'num_threads', 1, 'int');
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub duk{
	my $rH_cfg             = shift;
	my $log                = shift;
	my $ncontam            = shift;
	my $contam             = shift;
	my $db                 = shift;
	my $infile             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$contam, $ncontam, $log]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl'],
      ['duk', 'moduleVersion.duk']
    ]) . ' && ';
		$cmd .= ' memtime ';
		$cmd .= ' gunzip -c ' . $infile . ' |' ;	
		$cmd .= ' duk';
		$cmd .= ' -o ' . $log;
		$cmd .= ' -n ' . $ncontam;
		$cmd .= ' -m ' . $contam;
		$cmd .= ' -k '.LoadConfig::getParam($rH_cfg, 'duk', 'k', 1, 'int');
		$cmd .= ' -s '.LoadConfig::getParam($rH_cfg, 'duk', 's', 1, 'int');
		$cmd .= ' -c '.LoadConfig::getParam($rH_cfg, 'duk', 'c', 1, 'int');
		$cmd .= ' ' . $db;		
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub splitBarcodes{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $barcodes           = shift;
	my $outfile            = shift;
	my $log                = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' barcodes.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --barcodes ' . $barcodes;
		$cmd .= ' --outfile ' . $outfile;
		$cmd .= ' --num_threads '. LoadConfig::getParam($rH_cfg, 'barcodes', 'num_threads', 1, 'int');
		$cmd .= ' --log ' . $log;

		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub removeUnpairedReads{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfilePaired      = shift;
	my $unpairedR1         = shift;
	my $unpairedR2         = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfilePaired]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' removeUnpaired.pl';
		$cmd .= ' --infile '. $infile;
		$cmd .= ' --outfile_paired ' . $outfilePaired;
		$cmd .= ' --outfile_1 ' . $unpairedR1;
		$cmd .= ' --outfile_2 ' . $unpairedR2;
		$cmd .= ' --num_threads '.LoadConfig::getParam($rH_cfg, 'remove_unpaired', 'num_threads', 1, 'int');
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub splitPairs{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfileR1          = shift;
	my $outfileR2          = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfileR1, $outfileR2]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' splitPairs.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --outfile_1 ' . $outfileR1;
		$cmd .= ' --outfile_2 ' . $outfileR2;
		$cmd .= ' --num_threads '.LoadConfig::getParam($rH_cfg, 'split_pairs', 'num_threads', 1, 'int');
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub generateQscoreSheet{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $prefix             = shift;
	my $log                = shift;
	my $outfile            = shift;
	my $barcodes           = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['fastx', 'moduleVersion.fastx'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'qscoreSheets.pl ';
		$cmd .= ' --fastq ' . $infile;
		$cmd .= ' --tmp ' . LoadConfig::getParam($rH_cfg, 'default', 'tmpDir', 1, 'dirpath') ;
		$cmd .= ' --prefix ' . $prefix;
		$cmd .= ' --suffix suffix';
		$cmd .= ' --log ' . $log;
		$cmd .= ' --outfile ' . $outfile;
		$cmd .= ' --phred ' . LoadConfig::getParam($rH_cfg, 'default', 'qual', 1, 'int') ;
		$cmd .= ' --barcodes ' . $barcodes;
		$cmd .= ' --num_threads '. LoadConfig::getParam($rH_cfg, 'qscore_sheet', 'num_threads', 1, 'int');
		
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub generateQscoreGraphSingle{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $prefix             = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['R', 'moduleVersion.R'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' qscorePlots.pl';
		$cmd .= ' --infile_1 ' . $infile;
		$cmd .= ' --name ' . $prefix;  
		$cmd .= ' --pdf ' . $outfile; 
		$cmd .= ' --display 1';
		$cmd .= ' --single';

		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub generateQscoreGraphPaired{
	my $rH_cfg             = shift;
	my $infileR1           = shift;
	my $infileR2           = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infileR1, $infileR2] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['R', 'moduleVersion.R'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' qscorePlots.pl';
		$cmd .= ' --infile_1 ' . $infileR1;
		$cmd .= ' --infile_2 ' . $infileR2;
		$cmd .= ' --name qual_stats';  
		$cmd .= ' --pdf ' . $outfile; 
		$cmd .= ' --display 1';
		$cmd .= ' --paired';

		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub cutReads{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $begin              = shift;
	my $end                = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' &&';
		$cmd .=	' memtime';
		$cmd .= ' cutFastqSeq.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --begin ' . $begin;
		$cmd .= ' --end ' . $end;
		$cmd .= ' --outfile ' . $outfile;

		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub flash{
	my $rH_cfg             = shift;
	my $infileR1           = shift;
	my $infileR2           = shift;
	my $prefix             = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infileR1, $infileR2] , [$outdir.'/assembly_complete/ncontam_nphix_trimmed.extendedFrags.fastq']);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['flash', 'moduleVersion.flash'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' flash.pl';
		$cmd .= ' --infile_1 ' . $infileR1;
		$cmd .= ' --infile_2 ' . $infileR2;
		$cmd .= ' --prefix ' . $prefix;
		$cmd .= ' --outdir ' . $outdir;
		$cmd .= ' --n ' . LoadConfig::getParam($rH_cfg, 'flash', 'sampling', 1, 'int');
		$cmd .= ' --m ' . LoadConfig::getParam($rH_cfg, 'flash', 'minOverlap', 1, 'int');
		$cmd .= ' --M ' . LoadConfig::getParam($rH_cfg, 'flash', 'maxOverlap', 1, 'int');
		$cmd .= ' --x ' . LoadConfig::getParam($rH_cfg, 'flash', 'percentMismatch', 1, 'float');
		$cmd .= ' --p ' . LoadConfig::getParam($rH_cfg, 'flash', 'phred', 1, 'int');
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'flash', 'num_threads', 1, 'int');

		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub removePrimers{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $revPrimer          = shift;
	my $fwdPrimer          = shift;
	my $outfile            = shift;
	my $outfileFailed      = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'itagsQC.pl';
		$cmd .= ' --infile ' . $infile;
		if($revPrimer ne 'null'){
			$cmd .= ' --primer_3_prime ' . $revPrimer;
			$cmd .= ' --length_3_prime ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'length3Prime', 1, 'int');	
		}if($fwdPrimer ne 'null'){
			$cmd .= ' --primer_5_prime ' . $fwdPrimer;
			$cmd .= ' --length_5_prime ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'length5Prime', 1, 'int') ;	
		}
		#$cmd .= ' --qscore_1 ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'qscore1');
		#$cmd .= ' --qscore_2 ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'qscore2');
		$cmd .= ' --outfile ' . $outfile;
		$cmd .= ' --outfile_failed ' . $outfileFailed;
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'num_threads', 1, 'int');
		$cmd .= ' --qual ' . LoadConfig::getParam($rH_cfg, 'default', 'qual', 1, 'int');
		#$cmd .= ' --lq_threshold ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'lq_threshold');
		$cmd .= ' --primer_mismatch ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'primerMismatch', 1, 'int');
		#$cmd .= ' --min_length ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'minlength');
		#$cmd .= ' --N ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'N');
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub itagsQC{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $revPrimer          = shift;
	my $fwdPrimer          = shift;
	my $outfile            = shift;
	my $outfileFailed      = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile] , [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'itagsQC.pl';
		$cmd .= ' --infile ' . $infile;
		if($revPrimer ne 'null'){
			$cmd .= ' --primer_3_prime ' . $revPrimer;
			$cmd .= ' --length_3_prime ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'length3Prime', 1, 'int');	
		}if($fwdPrimer ne 'null'){
			$cmd .= ' --primer_5_prime ' . $fwdPrimer;
			$cmd .= ' --length_5_prime ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'length5Prime', 1, 'int') ;	
		}
		$cmd .= ' --qscore_1 ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'qscore1', 1, 'int');
		$cmd .= ' --qscore_2 ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'qscore2', 1, 'int');
		$cmd .= ' --outfile ' . $outfile;
		$cmd .= ' --outfile_failed ' . $outfileFailed;
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'num_threads', 1, 'int');
		$cmd .= ' --qual ' . LoadConfig::getParam($rH_cfg, 'default', 'qual', 1, 'int');
		$cmd .= ' --lq_threshold ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'lq_threshold', 1, 'int');
		$cmd .= ' --primer_mismatch ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'primerMismatch', 1, 'float');
		$cmd .= ' --min_length ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'minlength', 1, 'int');
		$cmd .= ' --N ' . LoadConfig::getParam($rH_cfg, 'itags_QC', 'N', 1, 'int');
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub countReport{
	my $rH_cfg             = shift;
	my $rA_files           = shift;
	my $rA_names           = shift;
	my $analysisType       = shift; 
	my $barcodesDist       = shift;
	my $OTUtable           = shift;
	my $obsTable           = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$OTUtable], [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'countReport.pl';
		foreach(@$rA_files){
			$cmd .= ' --file ' .$_;
		}
		foreach(@$rA_names){
			$cmd .= ' --name ' .$_;
		}
		$cmd .= ' --analysisType ' . $analysisType;
		$cmd .= ' --barcodesDist '. $barcodesDist;
		$cmd .= ' --OTUtable ' . $OTUtable;
		$cmd .= ' --obsTable ' . $obsTable;
		$cmd .= ' > ' .$outfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub txtToPdf{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $outfile            = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], [$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= 'txtToPdf.pl';
		$cmd .= ' --infile ' . $infile;
		$cmd .= ' --outfile ' . $outfile;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub mergePdf{
	my $rH_cfg             = shift;
	my $command            = shift;
	
	my $dummyOutfile       = "mergepdf.mugqic.done";
	my $ro_job = new Job();
	$ro_job->testInputOutputs([""] , [$dummyOutfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['ghostscript', 'moduleVersion.ghostscript'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' ' . $command;
		$cmd .= ' && touch ' . $dummyOutfile; 

		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub clustering1{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $barcodes           = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], ["$outdir/obs_filtered.fasta", "$outdir/obs_filtered.tsv"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['usearch', 'moduleVersion.usearch'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' clustering1.pl';
		$cmd .= ' --infile_fastq ' . $infile;
		$cmd .= ' --ref_db ' . LoadConfig::getParam($rH_cfg, 'DB', 'chimeras', 1, 'filepath'); 
		$cmd .= ' --barcodes ' . $barcodes;
		$cmd .= ' --outdir ' . $outdir;
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'clustering', 'num_threads', 1, 'int');
		#$cmd .= ' --start_at 4';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub clustering2{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $barcodes           = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], ["$outdir/obs_filtered.fasta", "$outdir/obs_filtered.tsv"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['usearch', 'moduleVersion.usearch'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' clustering2.pl';
		$cmd .= ' --infile_fastq ' . $infile;
		$cmd .= ' --ref_db ' . LoadConfig::getParam($rH_cfg, 'DB', 'chimeras', 1, 'path'); 
		$cmd .= ' --barcodes ' . $barcodes;
		$cmd .= ' --outdir ' . $outdir;
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'clustering', 'num_threads', 1, 'int');
		#$cmd .= ' --start_at 4';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub clustering3{
	my $rH_cfg             = shift;
	my $infile             = shift;
	my $barcodes           = shift;
	my $outdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile], ["$outdir/obs_filtered.fasta", "$outdir/obs_filtered.tsv"]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime'],
      ['usearch', 'moduleVersion.usearch'],
      ['dnaclust', 'moduleVersion.dnaclust'],
      ['tools', 'moduleVersion.tools'],
      ['perl', 'moduleVersion.perl']
    ]) . ' && ';
		$cmd .=	' memtime ';
		$cmd .= ' clustering3.pl';
		$cmd .= ' --infile_fastq ' . $infile;
		$cmd .= ' --ref_db ' . LoadConfig::getParam($rH_cfg, 'DB', 'chimeras', 1, 'filepath'); 
		$cmd .= ' --barcodes ' . $barcodes;
		$cmd .= ' --outdir ' . $outdir;
    $cmd .= ' --lowAbunCutOff ' . LoadConfig::getParam($rH_cfg, 'clustering', 'lowAbunCutOff', 1, 'int');
		$cmd .= ' --num_threads ' . LoadConfig::getParam($rH_cfg, 'clustering', 'num_threads', 1, 'int');
		#$cmd .= ' --start_at 4';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub clientReport{
  my $rH_cfg        = shift;
  my $iniFilePath   = shift;
  my $projectPath   = shift;
  my $pipelineType  = shift;
  my $reportPath    = shift;
    
  my $pipeline = 'pipeline=\"' .$pipelineType .'\",';
	my $titleTMP = LoadConfig::getParam($rH_cfg, 'report','projectName');
	my $title = 'report.title=\"' .$titleTMP .'\",';
	my $authorTMP = LoadConfig::getParam($rH_cfg, 'report','report.author');
	my $author = 'report.author=\"' .$authorTMP .'\",';
	my $contactTMP = LoadConfig::getParam($rH_cfg, 'report','report.contact');
	my $contact = 'report.contact=\"' .$contactTMP .'\",';

	my $ro_job = new Job();
	#$ro_job->testInputOutputs([$iniFilePath],[$projectPath]]);
	$ro_job->setUp2Date(0);

	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['R', 'moduleVersion.R']
    ]) . ' && ';
		$cmd .= ' R --no-save -e \'library(gqSeqUtils) ;';
		$cmd .= ' mugqicPipelineReport(';
		$cmd .= ' pipeline=\"' . $pipelineType . '\",';
		$cmd .= ' report.path=\"' . $reportPath . '\",';
		$cmd .= ' ini.file.path=\"' . $iniFilePath . '\",' ;
		$cmd .= ' report.title=\"' . LoadConfig::getParam($rH_cfg, 'report','projectName') . '\",' ;
		$cmd .= ' report.author=\"' . LoadConfig::getParam($rH_cfg, 'report','report.author') . '\",' ;
		$cmd .= ' report.contact=\"' . LoadConfig::getParam($rH_cfg, 'report','report.contact') . '\",' ;
		$cmd .= ' project.path=\"' . $projectPath . '\")\'' ;

		$ro_job->addCommand($cmd);
	}

	return $ro_job;
}

sub cleanup{
	my $rH_cfg             = shift;
	my $tmpdir             = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([""] , [""]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd = '';
    $cmd .= LoadConfig::moduleLoad($rH_cfg, [
      ['memtime', 'moduleVersion.memtime']
    ]) . ' && ';
		$cmd .= ' memtime ';
		$cmd .= 'rm '.$tmpdir.' -rf';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub templateSub{
	my $rH_cfg             = shift;
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
