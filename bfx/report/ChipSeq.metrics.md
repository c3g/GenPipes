### Sequencing and Alignment Metrics per Sample

General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: Sequencing and Alignment Statistics per Sample ([download table](trimMemSampleTable.tsv))

$trim_mem_sample_table$

* Raw Reads #: total number of reads obtained from the sequencer
* Surviving Reads #: number of remaining reads after the trimming step
* Surviving %: Surviving reads / Raw reads
* Aligned Filtered Reads: number of aligned reads to the reference after filtering by mapping quality
* Aligned Filtered %: Aligned Filtered Reads / Surviving Reads
* Duplicate Reads: number of aligned reads having the same 5' alignment positions (for both mates in the case of paired-end reads) after filtering by mapping quality
* Duplicate %: Duplicate / Aligned Filtered Reads
