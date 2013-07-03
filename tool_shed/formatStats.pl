#!/usr/bin/perl

# Outputs statistics for samples
# Louis Letourneau, Maxime Caron Jan 2012
# $ARGV[0] - Sample name regular expression

use strict;

my $EXPERIMENT_TYPE = 'wholegenome';

&main();

sub main {
  my @sampleNames = @ARGV;
  my @printedBasePctTitles;
  for my $sampleName (@sampleNames) {
	#print $sampleName . "\n";
    if(substr($sampleName,length($sampleName)-1,1) eq "/") {
      $sampleName = substr($sampleName,0,length($sampleName)-1);
    }

    my $directory = $sampleName."/";
    my $bamFile = $directory.$sampleName.'.sorted.dup.bam';
    my $insertSizeStatsFile;
    my $insertSizeStatsFile = $directory.$sampleName.'.sorted.dup.all.metrics.insert_size_metrics';
    my $flagStatsFile = $directory.$sampleName.'.sorted.dup.bam.flagstat';
    my $coverageFile = "";
    if($EXPERIMENT_TYPE eq "wholegenome") {
     $coverageFile = $directory.$sampleName.'.sorted.dup.all.coverage.sample_summary';
    }
    else {
     $coverageFile = $directory.$sampleName.'.sorted.dup.bam.CCDS.coverage.sample_summary'; 
    }

    if(! -e $insertSizeStatsFile) {
      print STDERR "$insertSizeStatsFile doesn't exist\n";
      next;
    }
    if(! -e $flagStatsFile) {
      print STDERR "$flagStatsFile doesn't exist\n";
      next;
    }
    if(! -e $coverageFile) {
      print STDERR "$coverageFile doesn't exist\n";
      next;
    }

    # Get lanes used in bam
    my $cmd = 'samtools view -H '.$bamFile.' | grep "^\@RG" | sed \'s/.*PU:run\([^_]\+_[0-9]\).*/\1/g\' | tr \'\n\' \' \' | sed \'s/ $//g\''; 
    #print STDERR "$cmd\n";
    my $lanesUsed = qx{$cmd};
    my @laneList = split(/ /, $lanesUsed);

    my $totalReads      = 0;
    my $nbQCPassedReads = 0;
    my $nbAlignedReads  = 0;
    my $uniqueReads     = 0;
    for my $lane (@laneList) {
      my $trimFile = '../reads/'.$sampleName.'/'.$lane.'/'.$sampleName.'.trim.out';
      open(COUNT, $trimFile) or die "Can't open $trimFile\n";
      while(my $line = <COUNT>) {
        if($line =~ /^Input Read Pairs/) {
          my ($inputReadNb, $survivingPairs) = $line =~ /^Input Read Pairs: (\d+) Both Surviving: (\d+) .*/;
          $totalReads += ($inputReadNb * 2);
          $nbQCPassedReads += ($survivingPairs * 2);
        }
      }
      close(COUNT);
    }

    my $nbDuplicateReads = 0;
    open(FLAGSTATS, "$flagStatsFile") or die "Can't open $flagStatsFile\n";
    while(my $line = <FLAGSTATS>) {
      chomp($line);
      if($line =~ /^([0-9]+) \+ [0-9]+ in total/) {
        if($1+0 != $nbQCPassedReads) {
          warn "Post trim read count doesn't match: ".$sampleName.': '.$1.' vs '.$nbQCPassedReads."\n";
          $nbQCPassedReads = $1+0;
        }
      }
      elsif($line =~ /^([0-9]+) \+ [0-9]+ duplicates/) {
        $nbDuplicateReads += $1;
      }
      elsif($line =~ /^([0-9]+) \+ [0-9]+ mapped/) {
        $nbAlignedReads += $1;
      }
    }
    close(FLAGSTATS);

    my $insertMedian = 0;
    my $insertMean = 0;
    my $insertMedianAbsDev = 0;
    my $insertStdDev = 0;
    open(INSERT, "$insertSizeStatsFile") or die "Can't open $insertSizeStatsFile\n";
    while(my $line = <INSERT>) {
      chomp($line);
      if($line =~ /^MEDIAN_INSERT_SIZE\tMEDIAN_ABSOLUTE_DEVIATION/) {
        $line = <INSERT>;
        chomp($line);
        my @values = split(/\t/, $line);
        $insertMedian = 0+$values[0];
        $insertMean = 0+$values[4];
        $insertMedianAbsDev = 0+$values[1];
        $insertStdDev = 0+$values[5];
        last;
      }
    }
    close(INSERT);

    my $meanCoverage;
    my @basePctTitles;
    my @basePctIdxs;
    open(COVERAGE, "$coverageFile") or die "Can't open $coverageFile\n";
    my $line = <COVERAGE>; # Header
    chomp($line);
    my @headerTitles = split(/\t/, $line);
    for (my $idx=0; $idx < @headerTitles; $idx++) {
      if($headerTitles[$idx] =~ /^%_bases_above_/) {
        push(@basePctTitles, $headerTitles[$idx]);
        push(@basePctIdxs, $idx);
      }
    }

    my $line = <COVERAGE>; # Sample Values
    chomp($line);
    my @values = split(/\t/, $line);
    $meanCoverage = 0+$values[2];
    my @pctCoverages;
    for (my $idx=0; $idx < @basePctTitles; $idx++) {
      push(@pctCoverages, $0+$values[ $basePctIdxs[$idx] ]);
    }
    close(COVERAGE);

    if(@printedBasePctTitles == 0) {
      @printedBasePctTitles = @basePctTitles;
      print "Sample Name,Lanes used (run_lane),Total,QC passed,Aligned,Uniquely Aligned,Non-Duplicates,Duplicate (%),Median Insert Size,Mean Insert Size,Average Deviation,Standard Deviation,Mean Coverage";
      for my $basePctTitle (@printedBasePctTitles) {
        print ",".$basePctTitle;
      }
      print "\n";
    }
    else {
      for (my $idx=0; $idx < @basePctTitles; $idx++) {
        if(!($basePctTitles[$idx] eq $printedBasePctTitles[$idx])) {
          die "Percent headers don't match between samples: ".join(',', @basePctTitles)." vs ".join(',', @printedBasePctTitles)."\n";
        }
      }
    }

    print $sampleName.",".$lanesUsed.",".$totalReads.",".$nbQCPassedReads.",".$nbAlignedReads.",NA,".($nbAlignedReads-$nbDuplicateReads).",".($nbDuplicateReads*100/$nbAlignedReads).",".$insertMedian.",".$insertMean.",".$insertMedianAbsDev.",".$insertStdDev.",".$meanCoverage;

    for my $pctCoverage (@pctCoverages) {
      print ",".$pctCoverage;
    }
    print "\n";
  }
}
