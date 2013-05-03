#!/bin/sh

##2012-03-12 - Mathieu Bourgey - mbourgey@genomequebec.com
# with help and use of scripts from Louis Letrouneau
##This pipeline allows random picking of 10000 R1 ends of PE reads and blastn them on the nt DB

#Usage runBlast.sh datafile($1) outputdir($2) paired($3)

#get args
FILES=(`cat $1`)
lengthF=`expr ${#FILES[@]} - 1`
OUTPUT_DIR=$2
PAIRED=$3

#params
FPR="/data/solexa/aligners/fastqPickRandom.pl --compressed --threshold"
FQ2FA="/data/solexa/aligners/fastq2FastaQual.pl"
BLAST="/data/solexa/aligners/ncbi-blast-2.2.25+/bin/blastn -query"
QSUB_OPTIONS="-cwd -V -N"

#check output directory
if [ ! -d $OUTPUT_DIR/ ]; then
        mkdir $OUTPUT_DIR
fi

if [ $PAIRED == 0 ]
then
	for i in `seq 0 1 $lengthF`
	do
		#get file path R1 
		j=`expr $i + 1`
		filepath=`echo ${FILES[$i]} | awk -F, '{print $1}'`
		sample=`echo ${FILES[$i]} | awk -F, '{print $2}'`
		fname=`basename $filepath`
		#Get number of sequences
		Nseq=$(zcat $filepath | awk ' { if ($0 == "+") { print $0} }' | wc -l) 
		echo "$sample has $Nseq sequences"
		#Get the threshold of random picking
		thrC=$(echo " scale=6; 10000 / $Nseq" | bc)
		echo "$sample has 0$thrC rdp threshold"
		#Random pick
		string="$FPR 0$thrC --input1 $filepath --out1 $OUTPUT_DIR/$sample.R1.RDP.fastq"
		echo $string
		echo $string | qsub -p 612 -j y -pe multicore 1 $QSUB_OPTIONS RDP_$sample 
		#Format fastq to fasta
		string="$FQ2FA $OUTPUT_DIR/$sample.R1.RDP.fastq $OUTPUT_DIR/$sample.R1.RDP.fasta $OUTPUT_DIR/$sample.R1.RDP.QUAL"
		echo $string
		echo $string | qsub -p 612 -j y -pe multicore 1 $QSUB_OPTIONS F2F_$sample -hold_jid RDP_$sample
		#Blast fasta
		string="$BLAST $OUTPUT_DIR/$sample.R1.RDP.fasta -db nt -out $OUTPUT_DIR/$sample.R1.RDP.blastres -perc_identity 80 -max_target_seqs 1 "
		echo $string
		echo $string | qsub -p 612 -j y -pe multicore 1 $QSUB_OPTIONS Blast_$sample -hold_jid F2F_$sample
		#subselect only the species and report only the 20 most frequent
		string="grep \">\" $OUTPUT_DIR/$sample.R1.RDP.blastres  | awk ' { print \$2 \"_\" \$3} ' | sort | uniq -c | sort -n -r | head -20 > $OUTPUT_DIR/$sample.R1.RDP.blastHit_20MF_species.txt"
		echo $string
		echo $string | qsub -p 612 -j y -pe multicore 1 $QSUB_OPTIONS fres_$sample -hold_jid Blast_$sample
	done
else
#run split and blast
	for i in `seq 0 2 $lengthF`
	do
		#get file path R1 and R2
		j=`expr $i + 1`
		filepath=`echo ${FILES[$i]} | awk -F, '{print $1}'`
		filepath2=`echo ${FILES[$j]} | awk -F, '{print $1}'`
		sample=`echo ${FILES[$i]} | awk -F, '{print $2}'`
		fname=`basename $filepath`
		#Get number of sequences
		Nseq=$(zcat $filepath | awk ' { if ($0 == "+") { print $0} }' | wc -l) 
		echo "$sample has $Nseq sequences"
		#Get the threshold of random picking
		thrC=$(echo " scale=6; 10000 / $Nseq" | bc)
		echo "$sample has 0$thrC rdp threshold"
		#Random pick
		string="$FPR 0$thrC --input1 $filepath --input2 $filepath2 --out1 $OUTPUT_DIR/$sample.R1.RDP.fastq --out2 $OUTPUT_DIR/$sample.R2.RDP.fastq"
		echo $string
		echo $string | qsub -p 612 -j y -pe multicore 1 $QSUB_OPTIONS RDP_$sample 
		#Format fastq to fasta
		string="$FQ2FA $OUTPUT_DIR/$sample.R1.RDP.fastq $OUTPUT_DIR/$sample.R1.RDP.fasta $OUTPUT_DIR/$sample.R1.RDP.QUAL"
		echo $string
		echo $string | qsub -p 612 -j y -pe multicore 1 $QSUB_OPTIONS F2F_$sample -hold_jid RDP_$sample
		#Blast fasta
		string="$BLAST $OUTPUT_DIR/$sample.R1.RDP.fasta -db nt -out $OUTPUT_DIR/$sample.R1.RDP.blastres -perc_identity 80 -max_target_seqs 1 "
		echo $string
		echo $string | qsub -p 612 -j y -pe multicore 1 $QSUB_OPTIONS Blast_$sample -hold_jid F2F_$sample
		#subselect only the species and report only the 20 most frequent
		string="grep \">\" $OUTPUT_DIR/$sample.R1.RDP.blastres  | awk ' { print \$2 \"_\" \$3} ' | sort | uniq -c | sort -n -r | head -20 > $OUTPUT_DIR/$sample.R1.RDP.blastHit_20MF_species.txt"
		echo $string
		echo $string | qsub -p 612 -j y -pe multicore 1 $QSUB_OPTIONS fres_$sample -hold_jid Blast_$sample
	done
fi