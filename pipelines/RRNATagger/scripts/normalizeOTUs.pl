#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Temp;
use Data::Dumper;
use Env qw/TMPDIR/;
use Statistics::R;
use File::Basename;
use Iterator::FastaDb;

$SIG{INT} = sub{exit}; #Handle ungraceful exits with CTRL-C.

my $usage=<<'ENDHERE';
NAME:
NormalizeOTUs.pl

PURPOSE:
To normalize a otu table. At the moment only quantile
normalization is implemented. TODO normalization to average.

INPUT:
--infile <string>         : File having cluster (OTUs) representative.
                            Headers must be properly formatted and have
                            barcodes information.

OUTPUT:
--outfile <string>        : outfile for OTU frequency table of all OTUs

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca
ENDHERE

## OPTIONS
my ($help, $infile, $infile_mapping, $outfile);
my $verbose = 0;

GetOptions(
	'infile=s'			=> \$infile,
	'outfile=s'			=> \$outfile,
    'help' 				=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile please provide at least one infile\n") 			unless($infile);
die("--outfile please provide an outfile for tabular values\n") unless($outfile);

## MAIN
no warnings 'numeric';
no warnings 'redefine';

my $tmpdir = File::Temp->newdir(
    "tmpDirItaggerNormalizeOTUsXXXXXX",
    DIR => $TMPDIR."/",
    CLEANUP => 0 #Let itagger.pl manage the removal of temp files.
);

## Use R to quantile normalize OTU table.
my $string = "";
my $R = Statistics::R->new();
$R->startR;
$R->run("options(stringsAsFactors = FALSE)");

$string .= '
## Sort and Rank each columns. For a matrix.
sortMatrix <- function(myCountsTable){
	myRankedTable = NULL
	for(i in 1:ncol(myCountsTable) ){
		myOriginalVector = as.numeric(myCountsTable[,i])
		mySortedVector = sort(as.numeric(myCountsTable[,i]))
		mySortedRankVector = NULL
		myUnsortedRankVector = NULL
		
		last = 999999999
		j = 999
		for(i in 1:length(mySortedVector)){
			if(i == 1){
				mySortedRankVector[i] = 1
				j = 1;
				last = mySortedVector[i]
			}else if(mySortedVector[i] == last){
				mySortedRankVector[i] = mySortedRankVector[i-1]
			}else{
				mySortedRankVector[i] = j + 1
				last = mySortedVector[i] 
				j = j + 1
			}
		}
	
		# Then the idea is that positions of originalVector and TempVector are identical.
		for(i in 1:length(myOriginalVector)){
			count = myOriginalVector[i]
			pos = match(count, mySortedVector)
			myUnsortedRankVector[i] = mySortedRankVector[pos]
			#mySortedVector[pos] = "NA"
		}
		#print("====")
		#print(myOriginalVector)
		#print(myUnsortedRankVector)
		#print("====")
		
		myRankedTable = cbind(myRankedTable, myUnsortedRankVector)
	}
	return (myRankedTable);
}

## Sort and Rank each columns. For a vector.
sortVector <- function(myCountsTable){
		myRankedTable = NULL
		myOriginalVector = as.numeric(myCountsTable)
		mySortedVector = sort(as.numeric(myCountsTable))
		mySortedRankVector = NULL
		myUnsortedRankVector = NULL
		
		#print("====")	
		#print(myOriginalVector)
		#print(myOriginalVector)
		#print(mySortedVector)
		#print("====")

		last = 999999999
		j = 999
		for(i in 1:length(mySortedVector)){
			if(i == 1){
				mySortedRankVector[i] = 1
				j = 1;
				last = mySortedVector[i]
			}else if(mySortedVector[i] == last){
				mySortedRankVector[i] = mySortedRankVector[i-1]
			}else{
				mySortedRankVector[i] = j + 1
				last = mySortedVector[i] 
				j = j + 1
			}
		}
	
		# Then the idea is that positions of originalVector and TempVector are identical.
		for(i in 1:length(myOriginalVector)){
			count = myOriginalVector[i]
			pos = match(count, mySortedVector)
			myUnsortedRankVector[i] = mySortedRankVector[pos]
			#mySortedVector[pos] = "NA"
		}
		print("====")
		print(myOriginalVector)
		print(myUnsortedRankVector)
		print("====")
		
		#myRankedTable = cbind(myRankedTable, myUnsortedRankVector)
		return (myUnsortedRankVector);
	
}

options(stringsAsFactors = FALSE)
data <- read.table("'.$infile.'", sep="\t", header=F, comment.char = "", skip=1)
myFirstLine = data[1,]
myOriginalHeader = data[1,]
myLineages = data[2:nrow(data),ncol(data)]
myClustersID = data[2:nrow(data),1]

myCountsTable = data[2:nrow(data),2:(ncol(data)-1)]
myRankedTable = sortMatrix(myCountsTable)
mySortedCounts = NULL

## Sort counts table
for(i in 1:ncol(myCountsTable) ){
	mySortedCounts = cbind(mySortedCounts, sort(as.numeric(myCountsTable[,i]), decreasing = FALSE) )
}
class(mySortedCounts) <- "numeric" ## Convert to numbers.

## Calculate row means and generate a new tables with means and their associated ranks.
myFinalRanks = NULL
myMeans = rowMeans(mySortedCounts, na.rm = FALSE, dims = 1)

## Remove any duplicates.
myMeans = unique(myMeans)
myMeansRank = sortVector(myMeans)
myFinalRanks = cbind(myFinalRanks, myMeans)
myFinalRanks = cbind(myFinalRanks, myMeansRank)

## To handle cases where value is 0...
max = max(myFinalRanks[,2])
for(i in 1:1000){
	myFinalRanks = rbind(myFinalRanks, c(0, (i+max) ) )
}

myNormalizedTable = NULL ## This is going to hold the final normalized table.
## Then generate a new table with the original myRankedTable
for(i in 1:ncol(myRankedTable)){
	myColValues = NULL;
	
	for(j in 1:nrow(myRankedTable) ){
		myValue = myRankedTable[j,i]
		myPos = match(myValue ,myFinalRanks[,2])

		#print( myValue )
		myMean = myFinalRanks[[myPos,1]]
		
		myColValues = append(myColValues, myMean)
	}
	myNormalizedTable = cbind(myNormalizedTable, myColValues)
}

## Add/restore headers to normalized table.
#myNormalizedTable = myNormalizedTable * 10000
#myNormalizedTable = round(myNormalizedTable, digits=0);
myNormalizedTable = cbind(myClustersID ,myNormalizedTable)
myNormalizedTable = cbind(myNormalizedTable, myLineages)
myNormalizedTable = data.frame(myNormalizedTable)
names(myNormalizedTable) = myOriginalHeader

## Outfiles to debug
#write.table(myFinalRanks, "C:/temp/myFinalRanks.tab", sep="\t")
#write.table(myRankedTable, "C:/temp/myRankedTable.tab", sep="\t")
#write.table(mySortedCounts, "C:/temp/mySortedCounts.tab", sep="\t")
write.table(myNormalizedTable, "'.$outfile.'.tmp", sep="\t", row.names = FALSE)

';

my @array;
push(@array, $string);

#$R->run('pdf("'.$outfile_graph.'" , width=7, height=7)');
#$R->run('par(mfrow=c('.$row.','.$col.'))'); 
foreach my $cmd (@array){
	$R->run($cmd);
	#print $cmd."\n";
}
#$R->run('dev.off()');

## FInally parse OTU table to remove " char and add a line to the header.
open(IN, "<".$outfile.".tmp") or die "Can't open file ".$outfile.".tmp\n";
open(OUT, ">".$outfile) or die "Can't open file ".$outfile."\n";
print OUT "#OTU TABLE - quantile normalized\n";
while(<IN>){
	chomp;
	my $line = $_;
	$line =~ s/"//g;
	print OUT $line."\n";
}
close(IN);
close(OUT);

## REMOVE TEMP FILES
sub END{
	system("rm ".$tmpdir." -rf");
}

1;
exit;
