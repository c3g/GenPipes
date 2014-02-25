#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp;
use Env qw/TMPDIR/;
use Statistics::R;
use File::Basename;

my $usage=<<'ENDHERE';
NAME:
PhylumBarplot.pl

PURPOSE:
Generates stacked barplots in pdf format
from multiple taxonomic summary tables 
typically generated with Qiime.

INPUT:
--infile <string>        : Can be multiple infiles.
--prefixes <string>         : Optional, if file names arent meaningful
                           , you can provide more meaningful names 
                           for each of the processed files. Important
                           , if this argument is provided, must be an
                           equal ammount of --infile and --prefixes args.
--same_color             : Set flag if triplicate in same color.
--bacteria_only          : set if graph is to be made with bacteria only.
--custom_ymax            : For developpement only, do not use.

OUTPUT:
--outfile_graph <string> : one pdf file contaning barplot graph.
--outfile_table <string> : one txt file containing tables.
NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, @infile, @prefixes, $outfile_graph, $outfile_table, $bacteria_only, $n, $title, $same_color, $custom_ymax);
my $verbose = 0;

GetOptions(
	'infile=s' 			=> \@infile,
	'prefixes=s' 		=> \@prefixes,
	'outfile_table=s' 	=> \$outfile_table,
	'outfile_graph=s' 	=> \$outfile_graph,
	'same_color'		=> \$same_color,
	'custom_ymax=i'		=> \$custom_ymax,
	'title=s'			=> \$title,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile please provide at least one infile\n") if(@infile < 1);
die("--outfile_table please provide an outfile for tabular values\n") unless($outfile_table);
die("--outfile_table please provide an outfile for pdf graph\n") unless($outfile_graph);
if(@prefixes){
	#die("Must be an equal ammount of --infile and --prefixes args") if(@prefixes != @infile);
}

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
	$col = 2;
	$cex = 0.7;
	$format = "c(1,10,20,30,40)";
}elsif($display == 2){
	$row = 4;
	$col = 2;
	$cex = 0.65;
	$format = "c(1,20,40)";
}

my @string;
#my @prefixes;
my $max=0;
my @max;
my $y_max=0;
my @xaxis = ();

my $R = Statistics::R->new();
$R->startR;

## First loop through all files to find maximum values of x and y.
my @max_file = @infile;
foreach my $max_file(@max_file){
	open(IN, $max_file) or die "Can't open file ".$max_file."\n";
	while(<IN>){
		chomp;
		next if($_ =~ m/\#/);
		if($_ =~ m/xmax/){
			my $curr_max = $_; 
			$curr_max =~ s/xmax: //;
			
			if($curr_max > $max){
				$max = $curr_max; 
				$R->run("xmax <- ".$max."");
			}
		}
		if($_ =~ m/series /){
			my $curr_max = 0; 
			$_ =~ s/series //;
			$_ =~ s/nan/0/g;
			my @row = split(/\t/, $_);
			foreach my $el (@row){
				if($el > $y_max){
					$y_max = $el; 
					$R->run("ymax <- ".$el."");
				}
			}
		}
	}
	close(IN);
}

## THEN BUILD GRAPH WITH R
my $curr_string = "";
$R->run("options(stringsAsFactors = FALSE)");
print $R->run("max") if($verbose);
$R->run("colors = c(
	'#191919',
	'#606060',
	'#A4A4A4',
	'#4C2A08',
	'#603912',
	'#996633',
	'#B36E35',
	'#C78F54',
	'#D3AD6C',
	'#E4C388',
	'#EFE9BD',
	'#6A1D1D',
	'#A40000',
	'#EB0000',
	'#FC4A2C',
	'#FF5E5E',
	'#191919',
	'#606060',
	'#A4A4A4',
	'#4C2A08',
	'#603912',
	'#996633',
	'#B36E35',
	'#C78F54',
	'#D3AD6C',
	'#E4C388',
	'#EFE9BD',
	'#6A1D1D',
	'#A40000',
	'#EB0000',
	'#FC4A2C',
	'#FF5E5E'
)");
#	'#FF8080',
#	'#FFA9A9',
#	'#EFBFBF',
#	'#A93617',
#	'#FF6000',
#	'#FE850F',
#	'#FEA70F',
#	'#FFE400',
#	'#FFF06C',
#	'#52007D',
#	'#A334AD',
#	'#C28CCF',
#	'#0E007A',
#	'#2A18B4',
#	'#645C99',
#	'#5C749C',
#	'#3878BA',
#	'#66A0DB',
#	'#8AC6E2',
#	'#103E0D',
#	'#056B00',
#	'#147B00',
#	'#30894E',
#	'#4BB65A',
#	'#65CD74',
#	'#50F324'

## ABSOLUTE COUNT
$curr_string .= '
	#my_mgp <- par()$mgp
	#my_mgp[1] <- 4
	#my_mgp[2] <- 0.5
	#my_mgp[3] <- 3
	#par(mgp=my_mgp)

	#my_mar <- par()$mar
	#my_mar[2] <- 4
	#par(mar=my_mar)
	#par(mfrow=c(2,2))

	#graph_data <- NULL
';

## DECLARE A NEW PLOT. WILL ADD POINTS IN NEXT STEPS.
#$curr_string .= '
#	plot(
#		c(0,xmax),
#		c(0,ymax),
#		pch=1,
#		cex=0.4,
#		lwd=0.1,
#		col=colors[0],
#		main="",
#		ylim=c(0,ymax),
#		xlim=c(0,xmax),
#		xlab="",
#		ylab="",
#		axes=F,
#	)
#';
open(OUT_TAB, ">".$outfile_graph.".tab") or die "Can't open file ".$outfile_graph."\n";

my $i=1;
my $k=0;
my $counter = 1;
my $init_i;
my @names;
# DRAW POINTS FOR EACH GRAPH.
foreach my $file (@infile){
	print $file."\n";
	open(IN, $file) or die "Can't open file ".$file."\n";
	my $prefix = shift(@prefixes) if(@prefixes);

	
	$init_i = $i;
	my @length = ();
	while(<IN>){
		chomp;
		next if($_ =~ m/\#/);
		
		if($_ =~ m/xaxis/){
			$_ =~ s/xaxis: //;
			my @row = split(/\t/, $_);
			$R->run("xaxis_".$init_i." <- NULL");
			
			my $length = @row;
			push(@length, $length);

			foreach my $value (@row){
     			$R->run("xaxis_".$init_i." <- c(xaxis_".$init_i.", ".$value.")");
			}

		}elsif($_ =~ m/series /){
			$_ =~ s/series //;
			$_ =~ s/nan/NA/g;
			$R->run("vector_".$i." <- NULL");
			my @row = split(/\t/, $_);
			my $length = @row;
			push(@length, $length);
			foreach my $value (@row){
				$R->run("vector_".$i." <- c(vector_".$i.", ".$value.")");
			}
			$i++;
		}elsif($_ =~ m/>> /){
			$_ =~ s/>> //;
			#if(@prefixes){
				push(@names, $prefix."-".$_);
			#}#else{
			#	push(@names, $_);
			#}
		}
	}
	close(IN);
	
	if($same_color){
		for(my $j=$init_i;$j<$i;$j++){
			$curr_string .= '
			lines(
				xaxis_'.$init_i.',
				vector_'.$j.',
				pch=20,
				lwd=0.5,
				cex=0.4,
				type="o",
				#ylim=c(0,10000),
				col=colors['.$counter.'],
				main="",
				xlab="",
				ylab="",
				#log="x",
				axes=F,
				#frame.plot=T,
				add=T
			)';
		}
	}else{
		my $k=0;
		for(my $j=$init_i;$j<$i;$j++){
			## Here we assume that the first 6 elements must be plotted on the same plot.	
			if($k == 0 || $k == 6){
				print "k:\t".$k."\n";	
				$curr_string .= '
					plot(
						c(0,xmax),
						c(0,ymax),
						pch=1,
						cex=0.4,
						lwd=0.1,
						col=colors[0],
						main="",
					';
				if($k == 0 && $custom_ymax){
					$curr_string .= '
						ylim=c(0,'.$custom_ymax.'),
					';
				}else{
					$curr_string .= '
						ylim=c(0,ymax),
					';
				}
				
				$curr_string .= '
						xlim=c(0,xmax),
						xlab="",
						ylab="",
						axes=F,
				)';
			}

			$curr_string .= '
			lines(
				xaxis_'.$init_i.',
				vector_'.$j.',
				pch=20,
				lwd=0.5,
				cex=0.4,
				type="o",
				col=colors['.($k+1).'],
				main="",
				xlab="",
				ylab="",
				axes=F,
				add=T
			)';
	
			if($k == 3 || $k == 6){
				$curr_string .= '
					box(lwd=0.8)
					axis(side=1, las=1, cex.axis=0.8,line=0, lwd=0.5)
					axis(side=2, las=1.5, cex.axis=0.8, line=0, lwd=0.5)
					mtext("Observed species",side=2,line=3, cex=1)
					mtext("Sequencing effort",side=1,line=3, cex=1)
					mtext("'.$prefix.'", side=3, line=1, cex = 0.6)
					abline(h=axTicks(side=2), col="lightgray", lty="dotted")		
					abline(v=axTicks(side=1), col="lightgray", lty="dotted")		
				';
					#abline(h='.$draw_line.', col="lightgray", lty=1)		
			}
			
			$k++;
		}	
	}

	$counter++;
}
close(OUT_TAB);


my $legend_string = "";

if($same_color){
	my $k=1;
	my $length = int(@names);
	foreach(@names){
		chomp;
		$legend_string .= '"'.$_.'"';
		$legend_string .= ',' if(($length - $k) != 0);
		$k++;
	}
}else{
	my $k=1;
	my $length = int(@names);
	foreach(@names){
		chomp;
		$legend_string .= '"'.$_.'"';
		$legend_string .= ',' if(($length - $k) != 0);
		$k++;
	}
}

print $legend_string."\n";

#	par(mfrow=c(1,1))
$curr_string .= '
	plot.new()
	legend(
		"bottom",
		legend = c('.$legend_string.'),
		fill = colors,
		ncol=2,
		cex=0.45,
		bty="n",
		border=T,
	)
';
my @array;
#print $curr_string."\n" if($verbose);
push(@array, $curr_string);
$R->run('pdf("'.$outfile_graph.'" , width=7, height=7)');
$R->run('par(mfrow=c('.$row.','.$col.'))'); 
foreach my $cmd (@array){
	$R->run($cmd);
	print $cmd."\n";
}
$R->run('dev.off()');

exit;
