#!/usr/bin/env perl

use warnings;
use strict;

use Env qw/PATH/;
use Getopt::Long;
use File::Copy;
use Cwd;
use Iterator::FastaDb;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use File::Which;

$| = 1;    # UNBUFFER OUTPUT

my $usage = <<'ENDHERE';
NAME:
Clustering.pl

PURPOSE:
To count occurrences of all sequences and optionally bin by abundance.
Intended for amplicon libraries. This script will cluster 16S rRNA
sequences in a fashion similar to R. Edgar's OTU pipe pipeline.
Briefly, there is: 1) 100% dereplication. 2) 99% clustering. 3)
removal of clusters < 3. 4) removal of chimeras. 5) Clustering at 97%
id.

INPUT:
--infile_fastq <fastq>     : Sequence file in fastq format.
--infile_fasta <fasta>     : Sequence file in fasta format.
--ref_db <string>          : Reference fasta database for chimera hunting.
--barcodes <string>        : Barcodes file in fasta format. Optional: used
                             for puting names in final cluster tables
                             instead of barcodes sequences.
--num_threads              : Number of threads. Only used in the 100%
                             dereplication step.

OUTPUT:
--outdir <string>          : Output directory where files will be generated.

OPTIONS:
--verbose                  : print status messages to stdout

AUTHOR/SUPPORT:
Julien Tremblay - julien.tremblay@mail.mcgill.ca

NOTE:
This clustering pipelne follows recommandations by Robert Edgar - drive5.com

ENDHERE

my ($help, $infile_fastq, $infile_fasta, $ref_db, $start_at, $verbose, $outdir, $barcodes, $num_threads);
GetOptions(
	'infile_fasta=s'	=> \$infile_fasta,
	'infile_fastq=s'	=> \$infile_fastq,
	'barcodes=s'		=> \$barcodes,
	'ref_db=s'			=> \$ref_db,
	'outdir=s'			=> \$outdir,
	'num_threads=i'		=> \$num_threads,
	'verbose'         	=> \$verbose,
	'start_at=i'		=> \$start_at,
	'help'            	=> \$help
);

if ($help){ print $usage; exit;}

## Validate
die "--infile_fastq or --infile_fasta arg missing.\n" if(!defined($infile_fastq) and !defined($infile_fasta));
die "--ref_db arg is missing (fasta reference database).\n" unless($ref_db);
die "--outdir arg missing.\n" unless($outdir);
$num_threads = 1 unless($num_threads);

## Check paths.
my $itaggerDereplicate = which('itaggerDereplicate.pl');
print "path of itaggerDereplicate:\t".$itaggerDereplicate."\n" if($verbose);
my $usearch = which('usearch');
print "path of usearch:\t".$usearch."\n" if($verbose);
$start_at = 0 unless($start_at);

##############
## BEG SUBS ##
##############
## Parse uclust table to generate a fasta file.
## Input uclust file
## Ref fasta file
## Fasta outfile
sub uclustToFasta{
	my ($uclust_file, $ref_fasta, $fasta_outfile) = @_;
	
	my %HoH = ();
	my %hash_seq = ();

	# Loop through uclust table.
	open(IN, $uclust_file) or die "Can't open ".$uclust_file."\n";
	while(<IN>){
		chomp;
		next if( substr($_,0,1) eq "#");
	
		my @row = split(/\t/, $_);
		my $status = $row[0];
		my $cluster_id = $row[1];
		my $query = $row[8];

		my @query = split(/;/, $query);
		my $seqid = shift(@query);
		my $size_field = pop(@query);
		$size_field =~ m/size=(\d+)/;			
		my $size = $1;

		#print "Size:\t".$size."\nStatus:\t".$status."\ncluster_ID:\t".$cluster_id."\nsize_field:\t".$size_field."\n" if($verbose);
	
		next if($status eq "C");	
		$HoH{$cluster_id}{size} += $size;
		$HoH{$cluster_id}{seqid} = $seqid if($status eq "S");
		
		# Parse barcodes
		foreach my $barcode_field (@query){
			$barcode_field =~ m/#([ACGT]*)/;
			my $barcode_seq = $1;
			$barcode_field =~ m/[ACGT]*=(\d+)/;
			my $barcode_count = $1;
			$HoH{$cluster_id}{$barcode_seq} += $barcode_count;
		}
	}
	close(IN);
	
	# Loop through fasta reference to grab actual sequence.
	my $ref_fasta_db = Iterator::FastaDb->new($ref_fasta) or die("Unable to open Fasta file ".$ref_fasta."\n");
	while( my $curr = $ref_fasta_db->next_seq() ) {
		my $header = $curr->header;
		$header =~ s/>//;
		my @row = split(/;/, $header);
		my $seqid = shift(@row);
		$hash_seq{$seqid} = $curr->seq;		
	}
	
	# Open fasta outfile for writing.
	open(OUT, ">".$fasta_outfile) or die "Can't open ".$fasta_outfile."\n";

	# Loop through hash and print.
	foreach my $cluster_id (sort {$HoH{$b}{size} <=> $HoH{$a}{size} } keys %HoH){
		print OUT ">".$cluster_id.";";
		
		for my $barcode_seq (keys %{ $HoH{$cluster_id} }){
			next if( ($barcode_seq eq "size") or ($barcode_seq eq "seqid") ); #if size or seqid
			print OUT "#".$barcode_seq."=".$HoH{$cluster_id}{$barcode_seq}.";";
		}
		print OUT "size=".$HoH{$cluster_id}{size}."\n";
		print OUT $hash_seq{ $HoH{$cluster_id}{seqid} }."\n";
	}
	close(OUT);	
}

## Write clusters table.
## Input: clusters fasta sequence rep.
## Output: cluster table outfile and cluster fasta file.
sub writeClusters{
	
	my ($clusters_fasta, $clusters_table, $new_cluster_fasta) = @_;

	my %HoH = ();
	my $i = 0;
	my %barcode_hash = ();
	
	## Loop to get all barcodes count.
	my $ref_fasta_db = Iterator::FastaDb->new($clusters_fasta) or die("Unable to open Fasta file ".$clusters_fasta);
	while( my $curr = $ref_fasta_db->next_seq() ) {
		my $header = $curr->header;
		my @row = split(/;/, $header);
		
		my $seqid = shift(@row);
		$seqid =~ s/>//;
		my $abundance = pop(@row);
	
		$abundance =~ m/(\d+)/;
		$HoH{$seqid}{abundance} = $1;
		
		foreach(@row){	
			if($_ =~ m/\#([ACGT]*)=(\d+)/){
				$HoH{$seqid}{$1} = $2;
				$barcode_hash{$1} = 1;
			}
		}
	}
	
	# Open outfile for writing.
	open(OUT, ">".$clusters_table) or die "Can't open ".$clusters_table."\n";
	
	# Print header
	my %names = ();
	if($barcodes){
		$ref_fasta_db = Iterator::FastaDb->new($barcodes) or die("Unable to open Fasta file ".$barcodes."\n");
		while( my $curr = $ref_fasta_db->next_seq() ) {
			my $header = $curr->header;
			$header =~ s/>//;
			$header =~ s/\t//g;
			$names{$curr->seq} = $header;
		}		
	}
	
	print OUT "#CLUSTER";
	foreach my $ref_barcode (sort keys %barcode_hash) {
		if($barcodes){
			#if(defined $names{$ref_barcode} and $names{$ref_barcode} ne "" and exists ($names{$ref_barcode}) ){
				print OUT "\t".$names{$ref_barcode};
			#}
		}else{
			print OUT "\t".$ref_barcode;
		}
	}
	print OUT "\n";
	
	# Then print table.
	foreach my $seqid (sort {$HoH{$b}{abundance} <=> $HoH{$a}{abundance} } keys %HoH){
		print OUT $seqid;	
	
		# Print barcodesin order.
		foreach my $ref_barcode (sort keys %barcode_hash) {
			if(exists $HoH{$seqid}{$ref_barcode} ){
				print OUT "\t".$HoH{$seqid}{$ref_barcode};
			}else{
				print OUT "\t0";
			}
		}
		print OUT "\n";        
	}
	close(OUT);

	# Print final clusters in seqobs format.	
	open(OUT, ">".$new_cluster_fasta) or die "Can't open ".$new_cluster_fasta." for writing\n";
	$ref_fasta_db = Iterator::FastaDb->new($clusters_fasta) or die("Unable to open Fasta file ".$clusters_fasta);
	while(my $curr = $ref_fasta_db->next_seq() ){
		$curr->header =~ m/(\d+)/;
		my $header = $1;
		print OUT ">".$header."\n".$curr->seq."\n";
	}
	close(OUT);

}
##############
## END SUBS ##
##############

## MAIN

## Dereplicate reads at 100%
my $cmd = $itaggerDereplicate;
$cmd .= " --fasta ".$infile_fasta if($infile_fasta);
$cmd .= " --fastq ".$infile_fastq if($infile_fastq);
$cmd .= " --outfile ".$outdir."/derep1.fasta";
$cmd .= " --num_threads ".$num_threads;
print $itaggerDereplicate."\n";
system($cmd) if($start_at <= 1);
die "command failed: $!\n" if($? != 0);

## Cluster at 99%
$cmd = $usearch." -cluster_smallmem ";
$cmd .= " ".$outdir."/derep1.fasta ";
$cmd .= " -id 0.99";
$cmd .= " -uc ".$outdir."/derep1_099.uc";
$cmd .= " -log ".$outdir."/derep1_099.log";
$cmd .= " -usersort";
system($cmd) if($start_at <= 2);
die "command failed: $!\n" if($? != 0);

## Parse Uclust file and write clusters
uclustToFasta($outdir."/derep1_099.uc", $outdir."/derep1.fasta", $outdir."/derep1_099.fasta") if($start_at <= 3);

## Sort by size
$cmd = $itaggerDereplicate;
$cmd .= " -fasta ".$outdir."/derep1_099.fasta";
$cmd .= " -outfile ".$outdir."/derep1_099_derep2.fasta";
$cmd .= " -minsize 3";
$cmd .= " -sort ";
system($cmd) if($start_at <= 3);
die "command failed: $!\n" if($? != 0);

#die $outdir."/derep1_099_derep2.fasta is empty.\n" if(-s $outdir."/derep1_099_derep2.fasta is empty.");

#Chimera detection (UCHIME de novo mode):
#  usearch -uchime amplicons.fasta [-chimeras ch.fasta] [-nonchimeras good.fasta]
#     [-uchimeout results.uch] [-uchimealns results.alns]
#  Input is estimated amplicons with integer abundances specified using ";size=N".
#usearch -uchime_denovo
$cmd = $usearch;
$cmd .= " -uchime_denovo ".$outdir."/derep1_099_derep2.fasta";
$cmd .= " -nonchimeras ".$outdir."/derep1_099_derep2_nonChimDeNovo.fasta";
#$cmd .= " -threads ".$num_threads;
system($cmd) if($start_at <= 4);
die "command failed: $!\n" if($? != 0);

#usearch -uchime_ref
#Chimera detection (UCHIME ref. db. mode):
#  usearch -uchime q.fasta [-db db.fasta] [-chimeras ch.fasta]
#    [-nonchimeras good.fasta] [-uchimeout results.uch] [-uchimealns results.alns]
$cmd = $usearch;
$cmd .= " -uchime_ref ".$outdir."/derep1_099_derep2_nonChimDeNovo.fasta";
$cmd .= " -db ".$ref_db;
$cmd .= " -strand both";
$cmd .= " -nonchimeras ".$outdir."/derep1_099_derep2_nonChimDeNovoRef.fasta";
$cmd .= " -threads ".$num_threads;
system($cmd) if($start_at <= 5);
die "command failed: $!\n" if($? != 0);

## Sort by size
$cmd = $itaggerDereplicate;
$cmd .= " -fasta ".$outdir."/derep1_099_derep2_nonChimDeNovoRef.fasta";
$cmd .= " -outfile ".$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted.fasta";
$cmd .= " -sort";
system($cmd) if($start_at <= 6);
die "command failed: $!\n" if($? != 0);

## Cluster at 97%
$cmd = $usearch." -cluster_smallmem";
$cmd .= " ".$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted.fasta";
$cmd .= " -id 0.97";
$cmd .= " -log ".$outdir."/derep1_099_derep2_nonChimDenovoRef_sorted_097.log";
$cmd .= " -uc ".$outdir."/derep1_099_derep2_nonChimDenovoRef_sorted_097.uc";
$cmd .= " -usersort";
system($cmd) if($start_at <= 7);
die "command failed: $!\n" if($? != 0);

## Parse Uclust file and write clusters
uclustToFasta(
	$outdir."/derep1_099_derep2_nonChimDenovoRef_sorted_097.uc",   # UC file 
	$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted.fasta",    # Ref file
	$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted_097.fasta" # new file
) if($start_at <= 8);

## Sort by size
$cmd = $itaggerDereplicate;
$cmd .= " --fasta ".$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted_097.fasta";
$cmd .= " --outfile ".$outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted_097_sorted.fasta";
$cmd .= " --sort";
system($cmd) if($start_at <= 9);
die "command failed: $!\n" if($? != 0);

## Generate Obs table
writeClusters($outdir."/derep1_099_derep2_nonChimDeNovoRef_sorted_097_sorted.fasta", $outdir."/obs_filtered.tsv", $outdir."/obs_filtered.fasta");

exit;
