
#!/bin/sh

# This script was written by Louis Letourneau and Maxime Caron
# Trims sequences for quality and length

if [ $FASTX_HOME ]; then
export FASTX_HOME=$FASTX_HOME
else
export FASTX_HOME=/data/solexa/tools/fastx_toolkit-0.0.13.1_modified/bin                   # Fastx home binary folder
fi

if [ $# != 5 ]; then
echo "usage: $0 file.fastq.gz quality length output_dir prefix"
exit
fi

echo "$1"
prefix=$5
zcat $1 | $FASTX_HOME/fastx_clipper -v -Q33 -l $3 -M11 -n -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | $FASTX_HOME/fastq_quality_trimmer -Q33 -t $2 -l $3 | gzip -c > "$4"${prefix}_t$2l$3.phred33.single.fastq.gz
