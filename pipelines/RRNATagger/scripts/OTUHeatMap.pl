#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp;
use Env qw/TMPDIR/;
use Statistics::R;
use File::Basename;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
OTUHeatMap.pl

PURPOSE:
Generates a heat map of OTUs. On the heat map, samples are 
displayed upon their distance as determined with hierarchical
clustering.

INPUT:
--infile <string>  : Can be multiple infiles (OTU table).
--n_most_abundant  : Only plot the n most abundant organims. TODO

OUTPUT:
--outdir <string> : one pdf file contaning barplot graph.
--prefix <string> : prefix to output file.

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $n, $outdir, $prefix);
my $verbose = 0;

GetOptions(
	'infile=s' 			    => \$infile,
	'outdir=s' 		      => \$outdir,
  'prefix=s'          => \$prefix,
	'n_most_abundant=i' => \$n,
  'help' 				      => \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile please provide at least one infile\n") unless($infile);
die("--outdir please provide an outdir for tabular values\n") unless($outdir);

my $legend_cex = 0.2;

## MAIN
no warnings 'numeric';
no warnings 'redefine';

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerOTUHeatMapXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

my $row = 1;
my $col = 1;
my $cex = 0.7;
my $format = "c(1,10,20,30,40)";

my $hist_string= "";
my @array;
my $counter = 0;
my $id;

my $tmp_infile = $tmpdir."/tmpinfile.txt";
my $number_of_samples = 0;

open(OUT, ">".$tmp_infile) or die "Can't open ".$tmp_infile."\n";
open(IN, "<".$infile) or die "Can't open file ".$infile."\n";
while(<IN>){
	chomp;

	# Get rid of the first line.
	if($counter == 0){
		$counter++;
		next;
	}
	
	# Calculate number of samples.
	if($counter == 1){
		my @row = split(/\t/, $_);
		$number_of_samples = @row;
	}
	
	# Then print all lines to temp file.
	my $line = $_;
	$line =~ m/#/;
	print OUT $line."\n";
	$counter++;
	
	if($n){
		last if($counter == $n);
	}
}
close(OUT);

## Set legend font size
if($counter < 10){
	$legend_cex = 0.40;
}elsif($counter < 20){
	$legend_cex = 0.40;
}elsif($counter < 40){
	$legend_cex = 0.30;
}elsif($counter < 80){
	$legend_cex = 0.2;
}elsif($counter < 160){
	$legend_cex = 0.15;
}

## Set graph font size
my $cex_names = 0.8;
if($number_of_samples < 10){
	$cex_names = 0.75;
}elsif($number_of_samples < 20){
	$cex_names = 0.65;
}elsif($number_of_samples < 40){
	$cex_names = 0.60;
}elsif($number_of_samples < 80){
	$cex_names = 0.55;
}elsif($number_of_samples < 97 ){
	$cex_names = 0.50;
}elsif($number_of_samples < 160){
	$cex_names = 0.4;
}

## BUILD GRAPH WITH R
my $R = Statistics::R->new();
$R->startR;
$R->run("options(stringsAsFactors = FALSE)");
$R->run("library(gplots)");

my $curr_string = '

	data <- read.table("'.$tmp_infile.'", sep="\t", comment.char="", header=T, skip=0, stringsAsFactors=F, row.names=1)
	mat=data.matrix(data[, 1:ncol(data)-1])
	coords <- which(mat == max(mat), arr.ind = TRUE)
	max <- mat[coords[,1], coords[,2]]
	
	names <- data[,-1]
	
#	palette_breaks <- lseq(1,max,length=100)
	palette_breaks <- exp(seq(log(1), log(max), length.out = 100))
	
	color_palette  <- colorRampPalette(c("#003366", "#6699FF", "#66CC00", "#FFFF33", "#FF0000"))(length(palette_breaks) - 1)
	heatmap.2(
		mat,
		Rowv=TRUE,
		Colv=TRUE,
		dendrogram= c("both"),
		density.info="histogram",
		distfun = dist,
		hclustfun = hclust,
		key=TRUE,
		keysize=1,
		trace="none",
		margins=c(10, 20),
		col=color_palette,
		labRow=data[,ncol(data)],
		breaks = palette_breaks,
		cexRow=0.4,
		cexCol='.$cex_names.'
	)
';

#$R->run("write.table(graph_data_perc, file='".$outfile_table.".perc', sep='\\t')");

push(@array, $curr_string);

## JPEG
$R->run('jpeg("'.$outdir.'/'.$prefix.'.jpeg", height=8, width=16, units="in", res=500)');
$R->run('par(mfrow=c('.$row.','.$col.'))'); 
foreach my $cmd (@array){
	$R->run($cmd);
	#print $cmd."\n";
}
$R->run('dev.off()');

## PDF
$R->run('pdf("'.$outdir.'/'.$prefix.'.pdf" , width=12, height=7)');
$R->run('par(mfrow=c('.$row.','.$col.'))'); 
foreach my $cmd (@array){
	$R->run($cmd);
	#print $cmd."\n";
}
$R->run('dev.off()');

## REMOVE TEMP FILES
sub END{
	system("rm ".$tmpdir." -rf");
}
1;
exit;
