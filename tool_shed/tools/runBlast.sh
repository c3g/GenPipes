#!/bin/sh

#get args
SAMPLE_TO=$1
OUTPUT_PREFIX=$2
FILE_1=$3
FILE_2=$4

OUTPUT_DIR=$(dirname ${OUTPUT_PREFIX})
PAIRED=0

if [ $# == 4 ]; then
  PAIRED=1
fi

#params
FPR="fastqPickRandom.pl --compressed --threshold"
FQ2FA="fastq2FastaQual.pl"
BLAST="blastn -query"

#check output directory
if [ ! -d $OUTPUT_DIR/ ]; then
  mkdir $OUTPUT_DIR
fi

Nseq=$(zcat $FILE_1 | awk ' { if ($0 == "+") { print $0} }' | wc -l) 
echo "$sample has $Nseq sequences"
#Get the threshold of random picking
thrC=$(echo " scale=6; $SAMPLE_TO / $Nseq" | bc)
echo "$sample has 0$thrC rdp threshold"

#Random pick
if [ $PAIRED == 0 ]
then
  $FPR 0$thrC --input1 $FILE_1 --out1 $OUTPUT_PREFIX.R1.subSampled_${SAMPLE_TO}.fastq
else
  $FPR 0$thrC --input1 $FILE_1 --input2 $FILE_2 --out1 $OUTPUT_PREFIX.R1.subSampled_${SAMPLE_TO}.fastq --out2 $OUTPUT_PREFIX.R2.subSampled_${SAMPLE_TO}.fastq
fi

#Format fastq to fasta
$FQ2FA $OUTPUT_PREFIX.R1.subSampled_${SAMPLE_TO}.fastq $OUTPUT_PREFIX.R1.subSampled_${SAMPLE_TO}.fasta $OUTPUT_PREFIX.R1.subSampled_${SAMPLE_TO}.qual

#Blast fasta
$BLAST $OUTPUT_PREFIX.R1.subSampled_${SAMPLE_TO}.fasta -db nt -out $OUTPUT_PREFIX.R1.subSampled_${SAMPLE_TO}.blastres -perc_identity 80 -num_descriptions 1 -num_alignments 1

#subselect only the species and report only the 20 most frequent
grep ">" $OUTPUT_PREFIX.R1.subSampled_${SAMPLE_TO}.blastres | awk ' { print $2 "_" $3} ' | sort | uniq -c | sort -n -r | head -20 > $OUTPUT_PREFIX.R1.RDP.blastHit_20MF_species.txt
