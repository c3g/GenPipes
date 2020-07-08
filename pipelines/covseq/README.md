Data Analysis Recommendations for MGI Sequencing Data
Paragon Genomics offers this panel of 343 amplicons covering 99.7% of the SARS-CoV-2 genome
(MN908947/NC_045512.2) with 92 bases uncovered at the ends of the genome. The 343 amplicons are
distributed into two pools. With amplicon size range from 116-196bp (median 149bp), the panel can be
used for DNBSEQ PE100 or PE150. Below are some notes on data processing.
We recommend the Broad Institute’s GATK Best Practice (https://software.broadinstitute.org/gatk/best-
practices/) as general guiding principles for sequencing data analysis.
1. Adapter Trimming.
With amplicon size ranges from 116 to 196bp, it is recommended to trim leftover adapter sequences for
PE150 sequencing (no need for PE100 sequencing) before read mapping.
Following is an example command with open source software cutadapt
(https://cutadapt.readthedocs.io/en/stable/).
cutadapt -g GAACGACATGGCTACGATCCGACTT \
-a AAGTCGGAGGCCAAGCGGTC \
-A AAGTCGGATCGTAGCCATGTCGTTC \
-G GACCGCTTGGCCTCCGACTT \
-e 0.1 -O 9 -m 20 -n 2 \
-o R1_out.fq.gz -p R2_out.fq.gz R1_in.fq.gz R2_in.fq.gz \
> cutadapt_report.output.txt
2. Map reads to reference genome.
The design was based on MN908947/NC_045512.2 and it is recommended to perform read mapping
against the reference sequence of NC_045512.2. Bwa mem is recommended for read mapping and de-
duplication procedure shall be skipped.
3. Trim primer sequences.
Before construction of a consensus genome sequence, it is recommended to remove primer sequences.
Software package fgbio is recommended. It requires primer genomic coordinates in a tab delimited file
which will be provided by Paragon Genomics to customers.
Following is an example command.
java -jar fgbio-1.2.0-e7ac607-SNAPSHOT.jar TrimPrimers -i input.bam -o
output.primerTrim.bam -p primer_info.tab -H true
UG4002-01
For Research Use Only. Not for use in diagnostic procedures.
41CleanPlex for MGI SARS-CoV-2 Research and Surveillance NGS Panel User Guide
4. Calculate QC metrics.
In order to assess the quality of the sequencing results, it is recommended to assign mapped reads to
amplicons based on mapping position. Subsequently, the following metrics can be used to measure
general performance of the panel.
● Mapping Rate: Percentage of reads mapped to reference genome. It assesses primer-dimers
and other PCR artifacts.
● On-Target Rate: Percentage of mapped reads that aligned to the targeted regions. It assesses
binding/amplification specificity of designed primers.
● Coverage Uniformity: Percentage of amplicons with read depth equal to or greater than 20% of
mean read depth of all amplicons in the panel. It measures performance uniformity of amplicons
in the panel.
To accommodate the calculation, a file in BED format listing amplicon start and end coordinates will be
provided to customers. The BED file can be downloaded from
www.paragongenomics.com/product_documents/.