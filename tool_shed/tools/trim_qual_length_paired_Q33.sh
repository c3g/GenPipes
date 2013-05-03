#!/bin/sh

# This script was written by Louis Letourneau and Maxime Caron
# Trims sequences for quality and length

if [ $FASTX_HOME ]; then
export FASTX_HOME=$FASTX_HOME
else
export FASTX_HOME=/lb/project/mugqic/software/tools/fastx_toolkit-0.0.13.1_modified/bin/                   # Fastx home binary folder
fi

if [ $PERL_FOLDER ]; then
export PERL_FOLDER=$PERL_FOLDER
else
export PERL_FOLDER=/lb/project/mugqic/software/bioinformatics/perl-tools/                         # Perl tool folder
fi


if [ $# != 6 ]; then
echo "usage: $0 file1.fastq.gz file2.fastq.gz quality length output_dir prefix"
exit
fi

prefix=$6

echo $1
echo $2


READ1_PRE=`mktemp -u /tmp/dmp-XXXXXX`;
mknod ${READ1_PRE} p || exit 1;
READ2_PRE=`mktemp -u /tmp/dmp-XXXXXX`;
mknod ${READ2_PRE} p || exit 1;
PAIR1_POST=`mktemp -u /tmp/dmp-XXXXXX`;
mknod ${PAIR1_POST} p || exit 1;
PAIR2_POST=`mktemp -u /tmp/dmp-XXXXXX`;
mknod ${PAIR2_POST} p || exit 1;
SINGLE_POST=`mktemp -u /tmp/dmp-XXXXXX`;
mknod ${SINGLE_POST} p || exit 1;

zcat $1 | $FASTX_HOME/fastx_clipper -v -Q33 -e -l 0 -M11 -n -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | $FASTX_HOME/fastq_quality_trimmer -Q33 -t $3 -l 0 > ${READ1_PRE} &
zcat $2 | $FASTX_HOME/fastx_clipper -v -Q33 -e -l 0 -M11 -n -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | $FASTX_HOME/fastq_quality_trimmer -Q33 -t $3 -l 0 > ${READ2_PRE} &
$PERL_FOLDER/split.pl ${READ1_PRE} ${READ2_PRE} ${PAIR1_POST} ${PAIR2_POST} ${SINGLE_POST} $4 &
cat ${PAIR1_POST}  | gzip -c > "$5"${prefix}_t$3l$4.phred33.pair1.fastq.gz &
cat ${PAIR2_POST}  | gzip -c > "$5"${prefix}_t$3l$4.phred33.pair2.fastq.gz &
cat ${SINGLE_POST} | gzip -c > "$5"${prefix}_t$3l$4.phred33.single.fastq.gz &
wait
rm ${READ1_PRE} ${READ2_PRE} ${PAIR1_POST} ${PAIR2_POST} ${SINGLE_POST}
