#!/usr/env/perl

=head1 NAME

I<ReadMetrics>

=head1 SYNOPSIS

ReadMetrics-> parseFlagstats()
ReadMetrics-> parseTrimOutput()
ReadMetrics-> mergeStats()

=head1 DESCRIPTION

B<ReadMetrics> is a library to generate QC, stats and metrics


=head1 AUTHOR
B<Johanna Sandoval> - I<johanna.sandoval@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package ReadMetrics;

# Strict Pragmas
#--------------------------
use strict;
use warnings;


#--------------------------

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------

sub parseFlagstats{
	my $sampleName    = shift;
	my $flagStatsFile = shift;
  my %stats         ;
  my $nbQCPassedReads = 0;
  my $nbDuplicateReads= 0;
  my $nbAlignedReads  = 0;
  open(FLAGSTATS, "$flagStatsFile") or die "Can't open $flagStatsFile\n";
  $stats{"_HEADER_"}=["Number of QC Passed Reads", "Number of Aligned Reads", "Number of Duplicate Reads"];
  while(my $line = <FLAGSTATS>) {
    chomp($line);
    if($line =~ /^([0-9]+) \+ [0-9]+ in total/) {
       $nbQCPassedReads = $1+0;
    }
    elsif($line =~ /^([0-9]+) \+ [0-9]+ duplicates/) {
      $nbDuplicateReads += $1;
    }
    elsif($line =~ /^([0-9]+) \+ [0-9]+ mapped/) {
      $nbAlignedReads += $1;
    }
  }
  $stats{$sampleName}=[$nbQCPassedReads, $nbAlignedReads, $nbDuplicateReads];
  close(FLAGSTATS);
  return \%stats;
}

sub parseTrimOutput{
  my $sampleName      = shift;
  my $trimOutputFile  = shift;
  my %stats           ;
  my $nbRawReads      = -1;
  my $nbFilteredReads = -1;
  my $nbSingleFiltered= -1;
  
  $stats{"_HEADER_"}=["Raw Fragments", "Fragment Surviving", "Single Surviving"];
  open(TRIMSTATS, "$trimOutputFile") or die "Can't open $trimOutputFile\n";
  while(my $line = <TRIMSTATS>) {
    chomp($line);
    if($line =~ /^Raw Fragments\,([0-9]+)/){
      $nbRawReads = $1;
    }elsif($line =~ /^Fragment Surviving\,([0-9]+)/ ){
      $nbFilteredReads = $1;
    }elsif($line =~ /^Single Surviving\,([0-9]+)/ ){
      $nbSingleFiltered= $1;
    }    
  }
  $stats{$sampleName}=[$nbRawReads, $nbFilteredReads, $nbSingleFiltered];
  close(TRIMSTATS);  
  return \%stats;
}

sub parseHomerAnnotations{
  my $annotationFile  = shift;
  my $outputFile      = shift;
  open(ANNOTATIONS, "$outputFile") or die "Can't open $outputFile\n";
  while(my $line = <TRIMSTATS>) {
    chomp($line);
    my @words=split( /\t/, $line);
    if($words =~ /^Raw Fragments\,([0-9]+)/){
      $nbRawReads = $1;
    }elsif($line =~ /^Fragment Surviving\,([0-9]+)/ ){
      $nbFilteredReads = $1;
    }elsif($line =~ /^Single Surviving\,([0-9]+)/ ){
      $nbSingleFiltered= $1;
    }    
  }

  

=for comment 
  
  a=`cat $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '{print $8}' | awk '$1 == "exon"' | wc -l`
  b=`cat $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '{print $8}' | awk '$1 == "intron"' | wc -l`
  c=`cat $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron"' | awk -F"\t" '$10>-2000 && $10<0' | wc -l`
  d=`cat $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron"' | awk -F"\t" '$10<=-2000 && $10>-10000' | wc -l`
  e=`cat $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron"' | awk -F"\t" '$10<=-10000 && $10>-100000' | wc -l`
  #f=`sed 1d $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron" && $8 !~ "TTS" && $8 !~ "promoter"' | awk -F'\t' '$10<-100000 || $10>100000' | wc -l`
  f=`sed 1d $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '$8 !~ "exon" && $8 !~ "intron"' | awk -F'\t' '$10<-100000 || $10>100000' | wc -l`
  #f=`cat $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '$10<-100000 || $10>100000' | wc -l`
  g=`sed 1d $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | wc -l`
  h=`expr $g - $a - $b - $c - $d - $e - $f`
  echo "exon,intron,proximal,distal,5d,gene_desert,other" > $OUTPUT_DIR/annotation/$designName/$designName.tss.stats.csv
  echo "$a,$b,$c,$d,$e,$f,$h" >> $OUTPUT_DIR/annotation/$designName/$designName.tss.stats.csv
  cat $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '{print $8}' | awk '$1 == "exon"' | awk -F',' '{print $2}' | sed -e 's/)//g' | awk '{print $2}' > $OUTPUT_DIR/annotation/$designName/$designName.exon.stats.csv
  cat $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '{print $8}' | awk '$1 == "intron"' | awk -F',' '{print $2}' | sed -e 's/)//g' | awk '{print $2}' > $OUTPUT_DIR/annotation/$designName/$designName.intron.stats.csv
  sed 1d $OUTPUT_DIR/annotation/$designName/$designName.annotated.csv | awk -F'\t' '{print $10}' > $OUTPUT_DIR/annotation/$designName/$designName.tss.distance.csv

=cut  
  
  
  $stats{"_HEADER_"}=["Raw Fragments", "Fragment Surviving", "Single Surviving"];
  open(TRIMSTATS, "$trimOutputFile") or die "Can't open $trimOutputFile\n";
  while(my $line = <TRIMSTATS>) {
    chomp($line);
    if($line =~ /^Raw Fragments\,([0-9]+)/){
      $nbRawReads = $1;
    }elsif($line =~ /^Fragment Surviving\,([0-9]+)/ ){
      $nbFilteredReads = $1;
    }elsif($line =~ /^Single Surviving\,([0-9]+)/ ){
      $nbSingleFiltered= $1;
    }    
  }
  $stats{$sampleName}=[$nbRawReads, $nbFilteredReads, $nbSingleFiltered];
  close(TRIMSTATS);  
  return \%stats;
}

sub mergeStats{
  my $sampleName           = shift;
  my $outputFile           = shift;
  my $rHoA_trimStats       = shift;
  my $rHoA_alignmentStats  = shift;
  my $printHeader          = shift;
  
  # Open an output file
  open(OUTFILE, ">".$outputFile);
  if (defined($printHeader) ){
    print OUTFILE $sampleName.",".join(',', @{$rHoA_trimStats->{"_HEADER_"}}). ','. join(',', @{$rHoA_alignmentStats->{"_HEADER_"}})."\n";
  }
	print OUTFILE $sampleName.",".join(',', @{$rHoA_trimStats->{$sampleName}}). ','. join(',', @{$rHoA_alignmentStats->{$sampleName}})."\n";
	close(OUTFILE);
}

1;