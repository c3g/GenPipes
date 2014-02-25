#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp;
use Env qw/TMPDIR/;
use Statistics::R;
use File::Basename;
use File::Which;

my $usage=<<'ENDHERE';
NAME:
itaggerPhylumBarplot.pl

PURPOSE:
Generates stacked barplots in pdf format
from multiple taxonomic summary tables 
typically generated with Qiime.

INPUT:
--infile <string>        : Can be multiple infiles.
--names <string>         : Optional, if file names arent meaningful
                           , you can provide more meaningful names 
                           for each of the processed files. Important
                           , if this argument is provided, must be an
                           equal ammount of --infile and --names args.
--bacteria_only          : set if graph is to be made with bacteria only.
--n_most_abundant        : Only plot the n most abundant organims.
--title <string>         : Graph title description. Optional.

OUTPUT:
--outfile_graph <string> : one pdf file contaning barplot graph.
--outfile_table <string> : one txt file containing tables.
NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, @infile, @names, $outfile_graph, $outfile_table, $bacteria_only, $n, $title);
my $verbose = 0;

GetOptions(
	'infile=s' 			=> \@infile,
	'names=s' 			=> \@names,
	'outfile_table=s' 	=> \$outfile_table,
	'outfile_graph=s' 	=> \$outfile_graph,
	'bacteria_only' 	=> \$bacteria_only,
	'n_most_abundant=i' => \$n,
	'title=s'			=> \$title,
    'help' 				=> \$help,
	'verbose'			=> \$verbose
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile please provide at least one infile\n") if(@infile < 1);
die("--outfile_table please provide an outfile for tabular values\n") unless($outfile_table);
die("--outfile_table please provide an outfile for pdf graph\n") unless($outfile_graph);
if(@names){
	#die("Must be an equal ammount of --infile and --names args") if(@names != @infile);
}

my $legend_cex = 0.2;
$title = " " unless($title);

## MAIN
no warnings 'numeric';
no warnings 'redefine';

my $display = 1;
my $cex;
my $format;
my $row;
my $col;
if($display == 1){
	$row = 2;
	$col = 1;
	$cex = 0.7;
	$format = "c(1,10,20,30,40)";
}elsif($display == 2){
	$row = 4;
	$col = 2;
	$cex = 0.65;
	$format = "c(1,20,40)";
}

my $tmpdir = $TMPDIR;
$tmpdir = $tmpdir."/" if(substr($tmpdir,-1) ne "/");
$tmpdir = $tmpdir."qscore_files/";
mkdir $TMPDIR unless -d $TMPDIR;
mkdir $tmpdir unless -d $tmpdir;

my $R_bin = which('R');
my $R = Statistics::R->new(r_bin=>$R_bin, log_dir=>$tmpdir, tmp_dir=>$tmpdir);
$R->startR;


my $hist_string= "";
my @array;
my $counter = 0;
my $id;

open(OUT_TAB, ">".$outfile_table);

## BUILD ONE TABLE FROM MULTIPLE L FILES
my %hash_all = ();
my @files = @infile;
my @array_of_hashref;
my @array_length;
my $header = "";
my $i=0;

foreach my $file (@files){
	open(IN, $file) or die "Can't open file ".$file."\n";
	while(my $line = <IN>){
		chomp($line);
		next if($. == 1);
		my @row = split(/\t/,$line);
		my $tax_id = $row[0];
		next if(!($tax_id =~ m/(k__bacteria|k__archaea)/i) and $bacteria_only);
		if(@row > 1){# In case something is wrong with the L tables...
			#print $tax_id."\n";
			$hash_all{$tax_id} = $tax_id;
		}else{
			next;
		}
	}
	close(IN);
}

foreach my $file (@files){	

	my %hash = ();
	open(IN, $file);
	my $entry_number = 0;

	while(<IN>){
		chomp($_);
		if($. == 1){
			my @header = split(/\t/, $_);
			shift(@header);
			$entry_number = scalar(@header);
			my $temp_header = join("\t", @header);
			$header .= $temp_header."\t";
			next;
		}
		
		my @row = split(/\t/,$_);
		my $length = scalar(@row);
		my $tax_id = $row[0];
		my $value = join("\t", @row[1..$length-1]);
	
		next if(!($tax_id =~ m/(k__bacteria|k__archaea)/i) and $bacteria_only);
		next if (@row <= 1); # In case something is wrong with the L tables.	
		if(exists $hash_all{$tax_id}){
			$hash{$tax_id} = $value."\t";
		}else{
			for(my $j=0; $j<$entry_number; $j++){
				$hash{$tax_id} = $hash{$tax_id}."0\t";
			}
		}
	}
	
	#fill with 0s absent taxonomies...
	while ( my ($key, $value) = each(%hash_all) ) {
		if(exists $hash{$key}){
	
		}else{
			$hash{$key} = "";
			for(my $j=0; $j<$entry_number; $j++){
				$hash{$key} = $hash{$key}."0\t";
			}
		}	
	}

	$array_of_hashref[$i] = \%hash;
	close(IN);
	$i++;
}
if(@names){
	my $named_header = "Lineage";
	for(my $k=0; $k<@names;$k++){
		$named_header .= "\t".$names[$k];
	}	
	print OUT_TAB $named_header."\n";
}else{
	print OUT_TAB "Lineage\t".$header."\n";
}

while ( my ($key, $value) = each(%hash_all) ) {
	print OUT_TAB $key."\t";
	$i=0;
	foreach(@array_of_hashref){
		if(exists $_->{$key}){ 
			#my $line = $_->{$key};
			#$line =~ s/\t\t/\t/gi
			print OUT_TAB $_->{$key};
		}else{ 
		
		}
		$i++;
	}
	print OUT_TAB "\n";
}
close(OUT_TAB);

## GET NAMES FOR GRAPHS
my $names_string = "";
my $number_of_samples = 0;

if(@names){
	foreach(@names){
		my $name = $_;
		$names_string .= "\'".$name."\',";
		$number_of_samples++;
	}
	chop($names_string);
}else{
	foreach my $file (@infile){
		#my $name = basename($_);
		open(FILE, $file) or die "Can't open file ".$file."\n";
		while(<FILE>){
			if($. == 1){
				my @row = split(/\t/, $_);
				shift(@row);
				foreach my $el (@row){
					push(@names, $el);
					$names_string .= "\'".$el."\',";
					$number_of_samples++;
				}
			}else{
				last;
			}
		}
		close(FILE);	
	}
	chop($names_string);
}
$names_string =~ s/\n//g;

my $names_table = '"Lineage",'.$names_string;

## Count number of different lineages:
my $lineages_number = 0;
my $command = "wc -l ".$outfile_table;
$lineages_number = `$command`;

## Set legend font size
if($lineages_number < 10){
	$legend_cex = 0.40;
}elsif($lineages_number < 20){
	$legend_cex = 0.40;
}elsif($lineages_number < 40){
	$legend_cex = 0.33;
}elsif($lineages_number < 80){
	$legend_cex = 0.33;
}elsif($lineages_number < 160){
	$legend_cex = 0.33;
}

## Set graph font size
my $cex_names = 1;
if($number_of_samples < 10){
	$cex_names = 0.8;
}elsif($number_of_samples < 20){
	$cex_names = 0.7;
}elsif($number_of_samples < 40){
	$cex_names = 0.6;
}elsif($number_of_samples < 80){
	$cex_names = 0.5;
}elsif($number_of_samples < 160){
	$cex_names = 0.3;
}elsif($number_of_samples < 200){
	$cex_names = 0.15;
}elsif($number_of_samples < 400){
	$cex_names = 0.05;
}

## BUILD GRAPH WITH R
$R->run("options(stringsAsFactors = FALSE)");
$R->run("data <- read.table('".$outfile_table."', sep='\t', header=F, skip=1)");
#$R->run("data <- data[order(data[,2]),]");
$R->run("colors = c(
	'#00FF00',
	'#FF8080',
	'#FF00FF',
	'#0000FF',
	'#00CCFF',
	'#CCFFFF',
	'#CCFFCC',
	'#99CCFF',
	'#CC99FF',
	'#FFCC99',
	'#3366FF',
	'#33CCCC',
	'#99CC00',
	'#FF99CC',
	'#FFCC00',
	'#FF9900',
	'#FF6600',
	'#666699',
	'#969696',
	'#003366',
	'#339966',
	'#003300',
	'#333300',
	'#FFFF99',
	'#993300',
	'#993366',
	'#333399',
	'#333333',
	'#000001',
	'#FFFFFF'
)");

## DEFINE APPROPRIATE LENGTH IT DEPENDS IF THERE ARE 
## MULTIPLE SAMPLES PER *L[\d+] TABLES.
my $matrix_length;
if(@names){
	$matrix_length = @names; 
}else{
	$matrix_length = @infile;
}

print STDERR "Starting to generate R graphs...\n" if($verbose);

## ABSOLUTE COUNT
my $curr_string = '
	pdf("'.$outfile_graph.'" , width=7, height=7);
	par(mfrow=c('.$row.','.$col.')); 
	#my_mgp <- par()$mgp
	#my_mgp[1] <- 4
	#my_mgp[2] <- 0.5
	#my_mgp[3] <- 3
	#par(mgp=my_mgp)

	my_mar <- par()$mar
	my_mar[2] <- 4
	par(mar=my_mar)
	#par(mfrow=c(2,2))

	names <- as.character(data[,1])
	graph_data <- NULL
	graph_data_perc <- NULL
	max_value <- 1;

	length <- '.@names.' + 1 
	for(i in 2:length){
		graph_data <- cbind(graph_data, as.numeric(data[,i]))
		curr_max_value <- sum(as.numeric(data[,i]))
		graph_data_perc <- cbind(graph_data_perc, (as.numeric(data[,i])/curr_max_value)*100)
		if( max_value <= curr_max_value ) max_value <- curr_max_value
	}
	graph_data_ref  <- graph_data
    
	if(ncol(graph_data_ref) == 1){
        graph_data      <- cbind(graph_data, graph_data_ref[,1:ncol(graph_data_ref)])
        graph_data_perc <- cbind(graph_data_perc, graph_data_ref[,1:ncol(graph_data_ref)])
    }else{
        graph_data      <- cbind(graph_data, apply(graph_data_ref[,1:ncol(graph_data_ref)],1,mean))
        graph_data_perc <- cbind(graph_data_perc, apply(graph_data_ref[,1:ncol(graph_data_ref)],1,mean))
    }

	row.names(graph_data)      <- names
	row.names(graph_data_perc) <- names

    if(ncol(graph_data_ref) == 1){
        graph_data      <- graph_data[order(-graph_data[,ncol(graph_data)]), ]
        graph_data_perc <- graph_data_perc[order(-graph_data_perc[,ncol(graph_data_perc)]), ]
    }else{
        graph_data      <- graph_data[order(-graph_data[,ncol(graph_data)]), ]
        graph_data_perc <- graph_data_perc[order(-graph_data_perc[,ncol(graph_data_perc)]), ]
	}
	
	
	names <- as.character(graph_data[,1])
	names_perc <- as.character(graph_data_perc[,1])

	graph_data <- as.matrix(graph_data)
	graph_data_perc <- as.matrix(graph_data_perc)
	
	';
	if($n){
		$curr_string .= '
		graph_data <- graph_data[1:'.$n.',]
		graph_data_perc <- graph_data_perc[1:'.$n.',]
		';
	}
	$curr_string .= '
	
	#----------------------------------------------
	op <- par(	
		oma = c(2,1,0,0) + 0.1,
		mar = c(2,2,1,1) + 0.1
	)

	barplot(
	    graph_data[, (1:ncol(graph_data)-1) ],
	    col=colors,
	    main="",
	    axes=F,
	    add=FALSE,
	    beside=FALSE,
	    las=2,
	    cex.names='.$cex_names.',
	    border=F,
	    ylim=c(0,max_value),
		names.arg=c('.$names_string.'),
	)
	box(lwd=0.4)
	axis(side=1, las=1, cex.axis=0.5, labels="", at=0, col="white", lwd=0.5)
	axis(side=2, las=1.5, cex.axis=0.5, line=0, lwd=0.5)
	mtext("Absolute abundance (number of reads)",side=2,line=3, cex=0.6)
	mtext("Taxonomy absolute '.$title.'" ,side=3,line=0.5, cex=0.6)
	#abline(h=axTicks(side=2), col="lightgray", lty="dotted")
	
	plot.new()
	legend(
		"bottom",
		legend = row.names(graph_data),
		fill = colors,
		cex='.$legend_cex.',
		ncol=4,
		bty="n",
		border=F,
	)

	barplot(
	     graph_data_perc[, (1:ncol(graph_data_perc)-1) ],
	     col=colors,
	     #log="y",
	     axes=F,
	     add=FALSE,
	     beside=FALSE,
	     #horiz=TRUE,
	     las=2,
	     cex.names='.$cex_names.',
	     border=F,
	     #yaxp=c(0, 100, 4),
	     ylim=c(0,100),
		 names.arg=c('.$names_string.'),
	     #legend.text=names,
	     #args.legend=list(x = "topright", border=F, bty="n", cex=0.9)
	)
	box(lwd=0.4)
	axis(side=1, las=1, cex.axis=0.5, labels="", at=0, col="white", lwd=0.5)
	axis(side=2, las=1.5, cex.axis=0.5, line=0, lwd=0.5)
	mtext("Relative abundance (%)",side=2,line=3, cex=0.6)
	mtext("Taxonomy relative '.$title.'" ,side=3,line=0.5, cex=0.6)
	#abline(h=axTicks(side=2), col="lightgray", lty="dotted")

	
	plot.new()
	legend(
		"bottom",
		legend = row.names(graph_data),
		fill = colors,
		cex='.$legend_cex.',
		ncol=4,
		bty="n",
		border=F,
	)
	write.table(graph_data, file="'.$outfile_table.'", sep="\\t", col.names=c('.$names_table.') )
	write.table(graph_data_perc, file="'.$outfile_table.'.perc", sep="\\t", col.names=c('.$names_table.') )

	#par(mfrow=c(1,1))
	#plot.new()
	#legend(
	#	"top",
	#	legend = row.names(graph_data),
	#	fill = colors,
	#	cex='.$legend_cex.',
	#	ncol=2,
	#	bty="n",
	#	border=F,
	#)
	dev.off();
';

$i++;
$counter++;
push(@array, $curr_string);

foreach my $cmd (@array){
	print STDERR $cmd."\n" if($verbose);
}

#$R->run('pdf("'.$outfile_graph.'" , width=7, height=7)');
print STDERR "Opened pdf...\n" if($verbose);
#$R->run('par(mfrow=c('.$row.','.$col.'))'); 
print STDERR "ran par command...\n" if($verbose);
foreach my $cmd (@array){
	print STDERR "Attempting to run R command...\n" if($verbose);
	$R->run($cmd);
	print STDERR "Ran command...\n" if($verbose);
	#print $cmd."\n" if($verbose);
}
#$R->run('dev.off()');

system "rm -rf ".$tmpdir; 
1;
exit;
