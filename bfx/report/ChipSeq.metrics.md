### Sequencing and Alignment Metrics per Sample

General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: Sequencing and Alignment Statistics per Sample ([download table](SampleMetrics.tsv))

$sample_table$

* Sample: Name of sample
* Mark Name: Name of Histone mark
* Raw Reads #: Total number of reads obtained from the sequencer
* Remaining Reads after Trimming #: Number of remaining reads after the trimming step
* Remaining Reads after Trimming %: 100 x Remaining Reads after Trimming / Raw reads
* Aligned Trimmed Reads #: Number of Aligned reads to the reference
* Aligned Trimmed Reads %: 100 x Aligned Trimmed Reads / Remaining Reads after Trimming
* Remaining Reads after Filtering #: Number of  remaining reads after filtering
* Remaining Reads after Filtering %: 100 x Remaining Reads after Filtering / Remaining Reads after Trimming
* Aligned Filtered Reads #: Number of Aligned reads to the reference after filtering
* Aligned Filtered Reads %: 100 x Aligned Filtered Reads / Filtered Remaining Reads
* Duplicate Reads #: Number of Duplicates reads (aligned reads having the same 5' alignment positions for both mates in the case of paired-end reads) after filtering
* Duplicate %: 100 x Duplicate Reads / Aligned Filtered Reads
* Final Aligned Reads # without Duplicates: Number of Aligned reads without duplicates count (Aligned Filtered Reads - Duplicate Reads)
* Mitchondrial Reads #: Number of reads aligned to chrM
* Mitochondrial %: 100 x Mitochondrial Reads / Aligned Filtered Reads 
