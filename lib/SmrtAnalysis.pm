#!/usr/env/perl

=head1 NAME

I<SmrtAnalysis>

=head1 SYNOPSIS

SmrtAnalysis->run()

=head1 DESCRIPTION

B<Picard> This a library to analyze PAcBio data using the SmrtAnalysis suite.

Input = file_name

Output = array

=head1 AUTHOR

Julien Tremblay - julien.tremblay@mail.mcgill.ca

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<Data::Dumper> Used to debbug

=cut

package SmrtAnalysis;

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
sub run {
	my $rH_cfg     = shift;
	my $paramsXml  = shift;
	my $inputXml   = shift;
	my $outdir     = shift;
	my $log        = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$paramsXml, $inputXml],[$log]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd;
		$cmd = 'module load ' . LoadConfig::getParam($rH_cfg, 'smrtanalysis','moduleVersion.smrtanalysis') .' &&'; 
		$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
		$cmd .= 'smrtpipe.py';
		$cmd .= ' -D NPROC=' . LoadConfig::getParam($rH_cfg, 'smrtanalysis', 'num_threads');
		$cmd .= ' -D TMP=' . LoadConfig::getParam($rH_cfg, 'smrtanalysis', 'tmpDir');
		$cmd .= ' --params=' .$paramsXml;
		$cmd .= ' --output=' .$outdir;
		$cmd .= ' --debug'; # ... Leave it there for now.
		$cmd .= ' xml:' .$inputXml; 
		$cmd .= ' > '. $log . ' 2>&1';
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub fofns {	
	my $rH_cfg     = shift;
	my $fofn       = shift;
	my $outputXml  = shift;
	my $outFofn    = shift;

  	my $ro_job = new Job();
	$ro_job->testInputOutputs([$fofn],[$outputXml]);

	if (!$ro_job->isUp2Date()) {
    	my $cmd;

		# Then do smrtpipe to generate input.xml
		$cmd = 'module load '. LoadConfig::getParam($rH_cfg, 'smrtanalysis', 'moduleVersion.smrtanalysis').' ;';
		$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
		$cmd .= 'fofnToSmrtpipeInput.py ';
		$cmd .= $fofn;
		$cmd .= ' > ' . $outputXml;
		$cmd .= ' && cp ' . $fofn . ' ' . $outFofn;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub referenceUploader{
 	
	my $rH_cfg     = shift;
	my $prefix     = shift;
	my $sampleName = shift;
	my $fasta      = shift;

  	my $ro_job = new Job();
	$ro_job->testInputOutputs([$fasta], );

  	if (!$ro_job->isUp2Date()) {
    	my $cmd;

	# Preload assembled contigs as reference
	$cmd = 'module load '. LoadConfig::getParam($rH_cfg, 'smrtanalysis', 'moduleVersion.smrtanalysis').' ;';
	$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
	$cmd .= ' referenceUploader';
	$cmd .= ' -c';
	$cmd .= ' -p ' . $prefix;
	$cmd .= ' -n ' . $sampleName;
	$cmd .= ' -f ' . $fasta;

    $ro_job->addCommand($cmd);
  }
  return $ro_job;
	
}

sub runCommand{
 	
	my $rH_cfg     = shift;
	my $string     = shift;	

  	my $ro_job = new Job();
	$ro_job->testInputOutputs(undef, undef );

  	if (!$ro_job->isUp2Date()) {
    	my $cmd;

		# Preload assembled contigs as reference
		$cmd = $string;
		#$cmd .= "cat ". LoadConfig::getParam($rH_cfg, 'default', 'polishingSettings');
	 	#$cmd .= " | sed \'s|REFERENCE|".$outdir."/".$sampleName.$suffix."/".$sampleName.$suffix." | g\' > ".$outdir."/".$sampleName.$suffix."/polishing/polishing.xml";


    $ro_job->addCommand($cmd);
  }
  return $ro_job;
	
}

sub blasr  {
	my $rH_cfg      = shift;
	my $infile      = shift;
	my $infileLong  = shift;
	my $outfile     = shift;
	my $outfileFofn = shift;
	my $sam         = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile, $infileLong],[$outfile, $outfileFofn]);
	
	if (!$ro_job->isUp2Date()) {
		#blasr filtered_subreads.fa filtered_longreads -out seeds.m4 -m 4 -nproc 8 -bestn 24 -nCandidates 24 -noSplitSubreads -minReadLength 200 -maxScore -1000 -maxLCPLength 16	
		my $cmd;
		$cmd = 'module load ' . LoadConfig::getParam($rH_cfg, 'default','moduleVersion.smrtanalysis') .' &&'; 
		$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
		$cmd .= 'blasr';
		$cmd .= ' ' . $infile;
		$cmd .= ' ' . $infileLong;
		$cmd .= ' -out ' . $outfile;
		$cmd .= ' -m ' . LoadConfig::getParam($rH_cfg, 'blasr','m');
		$cmd .= ' -nproc ' . LoadConfig::getParam($rH_cfg, 'blasr', 'num_threads');
		$cmd .= ' -bestn ' . LoadConfig::getParam($rH_cfg, 'blasr', 'bestn');
		$cmd .= ' -nCandidates ' . LoadConfig::getParam($rH_cfg, 'blasr', 'nCandidates');
		$cmd .= ' -noSplitSubreads'; 
		$cmd .= ' -minReadLength ' . LoadConfig::getParam($rH_cfg, 'blasr', 'minReadLength');
		$cmd .= ' -maxScore ' . LoadConfig::getParam($rH_cfg, 'blasr', 'maxScore');
		$cmd .= ' -maxLCPLength ' . LoadConfig::getParam($rH_cfg, 'blasr', 'maxLCPLength');
		if($sam eq "sam"){
			$cmd .= ' -sam';
		}
		$cmd .= ' && echo ' . $outfile . ' > ' .$outfileFofn;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub m4topre  {
	my $rH_cfg     = shift;
	my $infile     = shift;
	my $allm4      = shift;
	my $subreads   = shift;
	my $outfile    = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile],[$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd;
		$cmd = 'module load ' . LoadConfig::getParam($rH_cfg, 'smrtanalysis','moduleVersion.smrtanalysis') .' &&'; 
		$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
		$cmd .= 'm4topre.py';
		$cmd .= ' ' . $infile;
		$cmd .= ' ' . $allm4;
		$cmd .= ' ' . $subreads;
		$cmd .= ' ' . LoadConfig::getParam($rH_cfg, 'm4topre', 'bestn');
		$cmd .= ' > ' . $outfile;
		
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub pbdagcon  {
	my $rH_cfg       = shift;	
	my $infile       = shift;
	my $outfile      = shift;
	my $outfileFastq = shift;

	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile],[$outfile, $outfileFastq]);

	if (!$ro_job->isUp2Date()) {
		my $cmd;
		$cmd = 'module load ' . LoadConfig::getParam($rH_cfg, 'smrtanalysis','moduleVersion.smrtanalysis') .' &&'; 
		$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
		$cmd .= 'pbdagcon';
		$cmd .= ' -a ';
		$cmd .= ' -j ' . LoadConfig::getParam($rH_cfg, 'pbdagcon', 'num_threads');
		$cmd .= ' ' . $infile;
		$cmd .= ' > ' . $outfile;
		$cmd .= ' && ';
		$cmd .= ' cat ' . $outfile . ' |';
		$cmd .= ' awk \'{if($0~/>/){sub(/>/,"@",$0);print;}else{l=length($0);q="";while(l--){q=q "9"}printf("%s\n+\n%s\n",$0,q)}}\'';
		$cmd .= ' > ' . $outfileFastq;
	
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub quiver  {
	my $rH_cfg          = shift;
	my $infile          = shift;
	my $ref             = shift;
	my $outfileFasta    = shift;
	my $outfileVariants = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile],[$outfileFasta]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd;
		$cmd = 'module load ' . LoadConfig::getParam($rH_cfg, 'smrtanalysis','moduleVersion.smrtanalysis') .' &&'; 
		$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
		$cmd .= 'quiver';
		$cmd .= ' -j ' . LoadConfig::getParam($rH_cfg, 'quiver', 'num_threads');
		$cmd .= ' ' . $infile;
		$cmd .= ' -r ' . $ref;
		$cmd .= ' -o ' . $outfileVariants;
		$cmd .= ' -o ' . $outfileFasta;
		
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub variantcaller {
	my $rH_cfg          = shift;
	my $infile          = shift;
	my $ref             = shift;
	my $outfileFasta    = shift;
	my $outfileVariants = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile],[$outfileFasta]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd;
		$cmd = 'module load ' . LoadConfig::getParam($rH_cfg, 'smrtanalysis','moduleVersion.smrtanalysis') .' &&'; 
		$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
		$cmd .= 'variantCaller.py';
		$cmd .= ' -j ' . LoadConfig::getParam($rH_cfg, 'quiver', 'num_threads');
		$cmd .= ' --algorithm quiver';
		$cmd .= ' --referenceFilename ' . $ref;
		$cmd .= ' --parameters best';
		$cmd .= ' -o ' . $outfileVariants;
		$cmd .= ' -o ' . $outfileFasta;
		$cmd .= ' ' . $infile;
		
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

sub samtoh5  {
	my $rH_cfg          = shift;
	my $infile          = shift;
	my $ref             = shift;
	my $outfile         = shift;
	
	my $ro_job = new Job();
	$ro_job->testInputOutputs([$infile],[$outfile]);
	
	if (!$ro_job->isUp2Date()) {
		my $cmd;
		$cmd = 'module load ' . LoadConfig::getParam($rH_cfg, 'smrtanalysis','moduleVersion.smrtanalysis') .' &&'; 
		$cmd .= ' source \${SEYMOUR_HOME}/etc/setup.sh ; ';
		$cmd .= 'samtoh5';
		$cmd .= ' ' . $infile;
		$cmd .= ' ' . $ref;
		$cmd .= ' ' . $outfile;
		
		$ro_job->addCommand($cmd);
	}
	return $ro_job;
}

#sub generateFastq{
#	my $rH_cfg     = shift;
#	my $infile     = shift;
#	my $outfile    = shift;
#	
#	my $ro_job = new Job();
#	$ro_job->testInputOutputs([$infile],[$outfile]);
#	
#	if (!$ro_job->isUp2Date()) {
#		my $cmd;
#		$cmd .= ' cat ' . $infile . ' |';
#		$cmd .= ' awk \'{if($0~/>/){sub(/>/,"@",$0);print;}else{l=length($0);q="";while(l--){q=q "9"}printf("%s\n+\n%s\n",$0,q)}}\'';
#		$cmd .= ' > ' . $outfile;
#	
#		$ro_job->addCommand($cmd);
#	}
#	return $ro_job;
#}

1;

# Simple pbdagcon workflow script.  Written for the benefit of running via
# smrtpipe so I can communicate pipe errors to the task.  We're overcoming
# the limitation of smrtpipe forcing tasks to run serially, enabling a new
# level of pipelining that's extremely efficient in an imperfect world ...
# However, direct file I/O is faster by default.
#
#tmp=${tmp-"/tmp"}
#
#trap "rm -f $tmp/aln.$$.pre" EXIT SIGINT
#
## generate pre-alignments to a tmp directory
#m4topre.py $mym4 $allm4 $subreads ${bestn-24} > $tmp/aln.$$.pre
#
## pipe it to consensus and generate fasta
#pbdagcon -a -j ${nproc-15} $tmp/aln.$$.pre | tee ${fasta-"corrected.fa"} | \
## generate a fastq
#awk '{if($0~/>/){sub(/>/,"@",$0);print;}else{l=length($0);q="";while(l--){q=q "9"}printf("%s\n+\n%s\n",$0,q)}}' > ${fastq-"corrected.fq"}
#
#
## check the status of each piped command and exit non-zero if found
#for exitval in ${PIPESTATUS[*]}
#do
#    if [ $exitval -gt 0 ]
#    then
#        exit $exitval
#    fi
#done
#
#
#exit 0;
