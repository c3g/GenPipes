#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp;
use Env qw/TMPDIR/;
use File::Basename;
use Statistics::R;

my $usage=<<'ENDHERE';
NAME:
Procrustes.pl

PURPOSE:
Generates Procrustes tables using PCoA tables.

INPUT:
--infile <string>        : Can be multiple infiles.
--names <string>         : Optional, if file names arent meaningful
                           , you can provide more meaningful names 
                           for each of the processed files. Important
                           , if this argument is provided, must be an
                           equal ammount of --infile and --names args.

OUTPUT:
--outfile_graph <string> : one pdf file contaning barplot graph.
--outfile_table <string> : one txt file containing tables.

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, @infile, @names, $outfile_graph, $outfile_table, $outdir, $start_at);
my $verbose = 0;

GetOptions(
	'infile=s' 			=> \@infile,
	'names=s' 			=> \@names,
	'outdir=s'			=> \$outdir,
	'outfile_table=s' 	=> \$outfile_table,
	'outfile_graph=s' 	=> \$outfile_graph,
	'start_at=i'		=> \$start_at,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile please provide at least one infile\n") if(@infile < 1);
die("--outfile_table please provide an outfile for tabular values\n") unless($outfile_table);
die("--outfile_table please provide an outfile for pdf graph\n") unless($outfile_graph);
if(@names){
	#die("Must be an equal ammount of --infile and --names args") if(@names != @infile);
}

foreach(@infile){
	die "File is empty or does not exists <".$_.">. Mispelled the file?\n" if(!-s $_ || !-e $_); 
}

## Advanced usage: --start_at. If want to skip recalculation of coordinates, specify --start_at 2
$start_at = 1 unless($start_at);

mkdir $outdir unless -d $outdir;

my @MC;
my @M2;
for(my $i=0; $i<=@infile;$i++){
	$MC[$i] = "";
	$M2[$i] = "";
}

my @x = @infile;
pop(@x);
my @x_names = @names;
pop(@x_names);

my @y = @infile;
shift(@y);
my @y_names = @names;
shift(@y_names);

my %hash;

my $title_string = "";
my @plot_names;
my @plot_names_x;
my @plot_names_y;

open(OUT, ">".$outfile_table) or die "Can't open file ".$outfile_table."\n" if($start_at <= 1);
$title_string = $title_string."X/Y\t";
my $m = 1;
for(my $i=0; $i<@y;$i++){
	my $file_y = $y[$i];
	my $base_y = basename($file_y);
	$base_y =~ s/\.tab//;
	$base_y =~ s/\.txt//;
	
	my $curr_y_name = $y_names[$i];

	for(my $j=0; $j<$m;$j++){

	my $curr_x_name = $x_names[$j];

		if($j == 0){
			$MC[$i+1] = $y_names[$m-1]."\t";
			$M2[$i+1] = $y_names[$m-1]."\t";
			#next;
		}
		
		my $file_x = $x[$j];
		my $base_x = basename($file_x);
		$base_x =~ s/\.tab//;
		$base_x =~ s/\.txt//;	
		
		#print $file_x."\n";

		# Run cmd
		my $cmd = "transform_coordinate_matrices.py -i ".$file_y.",".$file_x." -r 10000 -o ".$outdir."/";
		system($cmd) if($start_at <= 1);
		print $cmd."\n";
		
		## Extract values from file.
		## FP1 FP2 Included_dimensions MC_p_value Count_better M^2
		open(PRO, "<".$outdir."/".$base_y."_".$base_x."_procrustes_results.txt") or die "Can't open file ".$outdir."/".$base_y."_".$base_x."_procrustes_results.txt";
		while(<PRO>){
			chomp;
			next if($. == 1);
			my @row = split(/\s/, $_);
			$MC[$i+1] = $MC[$i+1]."".$row[3]."\t";
			$M2[$i+1] = $M2[$i+1]."".$row[5]."\t";
		}

		## rename pc1_transformed.txt and pc2_transformed.txt
		system("mv ".$outdir."/pc1_transformed.txt ".$outdir."/pc1_transformed_".$curr_y_name.".txt") if($start_at <= 1);
		system("mv ".$outdir."/pc2_transformed.txt ".$outdir."/pc2_transformed_".$curr_x_name.".txt") if($start_at <= 1);
		$hash{$curr_y_name."+".$curr_x_name}{pc1} = $outdir."/pc1_transformed_".$curr_y_name.".txt";
		$hash{$curr_y_name."+".$curr_x_name}{pc2} = $outdir."/pc2_transformed_".$curr_x_name.".txt";
		$hash{$curr_y_name."+".$curr_x_name}{pc1_name} = $curr_y_name;
		$hash{$curr_y_name."+".$curr_x_name}{pc2_name} = $curr_x_name;
		#push(@plot_names, $curr_y_name."+".$curr_x_name);
		close(PRO);	
	}	
	$m++;
}

if($start_at <= 1){
	## Print M2 and Monte-Carlo values.
	print OUT "\nMonte-Carlo\n";
	print OUT $title_string;
	for(@MC){
		print OUT $_."\n";
	
	}
	
	print OUT "\nM^2\n";
	print OUT $title_string;
	foreach(@M2){
		print OUT $_."\n";
	
	}
}

## PRINT PLOTS with R.
my $R = Statistics::R->new();
$R->startR;
my $curr_string = '';
$R->run("options(stringsAsFactors = FALSE)");
#$R->run("colors = c(
#	'#990000',
#	'#FF0000',
#	'#FF6666',
#	'#009900',
#	'#00FF00',
#	'#99FF99',
#	'#000099',
#	'#0066CC',
#	'#3399FF',
#	'#990099',
#	'#FF3399',
#	'#FF66B2',
#	'#000000',
#	'#202020',
#	'#404040',
#	'#FF2400',
#	'#FF7f00',
#	'#FF4500',
#	'#00FF00',
#	'#FF8080',
#	'#FF00FF',
#	'#0000FF',
#	'#00CCFF',
#	'#CCFFFF',
#	'#CCFFCC',
#	'#99CCFF',
#	'#CC99FF',
#	'#FFCC99',
#	'#3366FF',
#	'#33CCCC',
#	'#99CC00',
#	'#FF99CC',
#	'#FFCC00',
#	'#FF9900',
#	'#FF6600',
#	'#666699',
#	'#969696',
#	'#003366',
#	'#339966',
#	'#003300',
#	'#333300',
#	'#FFFF99',
#	'#993300',
#	'#993366',
#	'#333399',
#	'#333333',
#	'#000001',
#	'#FFFFFF'
#)");

#$R->run("colors = c(
#	'#990000',
#	'#FF0000',
#	'#FF6666',
#	'#009900',
#	'#00FF00',
#	'#99FF99',
#	'#000099',
#	'#0066CC',
#	'#3399FF',
#	'#990099',
#	'#990000',
#	'#FF0000',
#	'#FF6666',
#	'#009900',
#	'#00FF00',
#	'#99FF99',
#	'#000099',
#	'#0066CC',
#	'#3399FF',
#	'#990099'
#	)"
#);


$R->run("colors = c(
	'#990000',
	'#009900',
	'#000099',
	'#990099',
	'#000000',
	'#FF2400',
	'#00FF00',
	'#0000FF',
	'#CCFFCC',
	'#FFCC99',
	'#990000',
	'#009900',
	'#000099',
	'#990099',
	'#000000',
	'#FF2400',
	'#00FF00',
	'#0000FF',
	'#CCFFCC',
	'#FFCC99'
)");
#	'#99CC00',
#	'#FF9900',
#	'#969696',
#	'#003300',
#	'#993300',
#	'#333333'

## Tweak margin appearance here.
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
';

my $counter = 0;
my $i = 1;
my $j = 0;
my $pch = 16;
while ( my ($key, $value) = each(%hash) ) {
		
	my @legend_names;
	my $sample_number = 0;
	my $pc1_name = $hash{$key}{pc1_name};
	my $pc2_name = $hash{$key}{pc2_name};
	my $plot_name = $key;
	my $color_counter = 0;
	my $color_counter_2 = 0;

	## There will be two files for each $pc keys.
	for my $pc (sort  keys %{ $hash{$key} } ) {
		next if($pc eq "pc1_name");
		next if($pc eq "pc2_name");

		print $key."\t".$hash{$key}{$pc}."\n";
		my $file = $hash{$key}{$pc}; 
		
		my @length = ();
		$R->run("vectorX_".$i." <- NULL");
		$R->run("vectorY_".$i." <- NULL");
		my $x_label = "";
		my $y_label = "";	
	
		open(IN, "<".$file) or die "Can't open file ".$file."\n";		
		while(<IN>){
			chomp;
			next if($_ =~ m/\#/);
			next if($_ =~ m/eigvals/);
			next unless length;
			next if($_ =~ m/pc vector number/);
		
			if($_ =~ m/% variation explained/){
				$_ =~ s/% variation explained\t//;
				my @row = split(/\t/, $_);
				$x_label = sprintf('%.2f', $row[0]);
				$y_label = sprintf('%.2f', $row[1]);
				$x_label = "P1 ".$x_label."%";
				$y_label = "P2 ".$y_label."%";
		
			}else{
				my @row = split(/\t/, $_);
				$R->run("vectorX_".$i." <- c(vectorX_".$i.", ".$row[1].")");
				$R->run("vectorY_".$i." <- c(vectorY_".$i.", ".$row[2].")");
				push(@legend_names, ($row[0]."_".$pc1_name)) if($j == 0);
				push(@legend_names, ($row[0]."_".$pc2_name)) if($j == 1);
				$sample_number++;
				$color_counter++;
			}	
		}
		close(IN);
		
		#$R->run("write.table(vectorX_".$i.", file='".$outdir."/vectorX_".$i.".txt')");
		#$R->run("write.table(vectorY_".$i.", file='".$outdir."/vectorY_".$i.".txt')");

		$pch=16 if($j == 0);
		$pch=7 if($j > 0);
		
		if($j == 0){

			$curr_string .= '
				plot(
					vectorX_'.$i.',
					vectorY_'.$i.',
					pch='.$pch.',
					cex=1.3,
					lwd=0.1,
					col=colors[1:'.$color_counter.'],
					main="'.$plot_name.'",
					cex.main=0.7,
					ylim=c(-0.35, 0.35),
					xlim=c(-0.35, 0.35),
					xlab="",
					ylab="",
					axes=F,
				)
				#text(vectorX_'.$i.', vectorY_'.$i.', paste(vectorX_'.$i.', vectorY_'.$i.', sep=", "), cex=0.3);
				box(lwd=0.8)
				axis(side=1, las=1, cex.axis=0.8,line=0, lwd=0.5)
				axis(side=2, las=1.5, cex.axis=0.8, line=0, lwd=0.5)
				mtext("'.$y_label.'",side=2,line=3, cex=1)
				mtext("'.$x_label.'",side=1,line=3, cex=1)
				abline(h=axTicks(side=2), col="lightgray", lty="dotted")		
				abline(v=axTicks(side=1), col="lightgray", lty="dotted")		
			';
			$color_counter_2 = $color_counter + 1;
			$j++;
		}else{ # Else means this is the second pc file.
			$pch = 6;
			$curr_string .= '
				#par(new=TRUE)
				points(
					vectorX_'.$i.',
					vectorY_'.$i.',
					pch='.$pch.',
					cex=1.3,
					lwd=2,
					col=colors['.($color_counter_2).':'.($color_counter).'],
					xlab="",
					ylab="",
					axes=F,
				)
				#text(vectorX_'.$i.', vectorY_'.$i.', paste(vectorX_'.$i.', vectorY_'.$i.', sep=", "), cex=0.3);
				#box(lwd=0.8)
				#axis(side=1, las=1, cex.axis=0.8,line=0, lwd=0.5)
				#axis(side=2, las=1.5, cex.axis=0.8, line=0, lwd=0.5)
				#mtext("'.$y_label.'",side=2,line=3, cex=1)
				#mtext("'.$x_label.'",side=1,line=3, cex=1)
				#abline(h=axTicks(side=2), col="lightgray", lty="dotted")		
				#abline(v=axTicks(side=1), col="lightgray", lty="dotted")		
			';	

			# Legend string is drawn here because we assume there is two points set per plot.	
			my $legend_string = "";
		
			my $k=1;
			my $length = int(@legend_names);
			foreach(@legend_names){
				chomp;
				$legend_string .= '"'.$_.'"';
				$legend_string .= ',' if(($length - $k) != 0);
				$k++;
			}
		
			$R->run("legend_".$i." <- c(".$legend_string.")");
		
			#par(mfrow=c(1,1))
			$curr_string .= '
				plot.new()
				legend(
					"bottom",
					legend = c(legend_'.$i.'),
					fill = colors,
					cex=0.3,
					ncol=2,
					bty="n",
					border=T,
				)
			';
			$j = 0;
			$color_counter = 0;
		}
		$i++;
	}
}

my @array;
my $row = 2;
my $col = 2;
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
