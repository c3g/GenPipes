
#!/bin/sh

# This script was written by Louis Letourneau and Maxime Caron
# Trims sequences for quality and length

export FASTX_HOME=/data/solexa/tools/fastx_toolkit-0.0.13.1_modified/bin          # Fastx home binary folder
export EMBOSS_HOME=/data/solexa/aligners/EMBOSS-6.2.0/bin                         # Emboss home binary folder

if [ $# != 5 ]; then
echo "usage: $0 file.fastq.gz quality length output_dir prefix"
exit
fi

echo "$1"
prefix=$5
zcat $1 | $FASTX_HOME/fastx_clipper -v -l $3 -M11 -n -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | $FASTX_HOME/fastq_quality_trimmer -t $2 -l $3 | $EMBOSS_HOME/seqret fastq-illumina::stdin fastq::stdout | gzip -c > "$4"${prefix}_t$2l$3.phred33.single.fastq.gz
