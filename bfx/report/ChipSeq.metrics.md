### Sequencing and Alignment Metrics per Sample

General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: Sequencing and Alignment Statistics per Sample ([download table](trimMemSampleTable.tsv))

$sample_table$

* Sample: Name of sample
* Mark Name: Name of Histone mark
* Raw Reads #: Total number of reads obtained from the sequencer
* Trimmed Remaining Reads #: Number of remaining reads after the trimming step
* Trimmed Remaining Reads %: 100 x Trimmed Remaining Reads / Raw reads
* Aligned Reads #: Number of Aligned reads to the reference
* Aligned Reads %: 100 x Aligned Reads / Trimmed Remaining Reads
* Filtered Remaining Reads #: Number of  remaining reads after filtering
* Filtered Remaining Reads %: 100 x Filtered Remaining Reads / Trimmed Remaining Reads
* Aligned Filtered Reads #: Number of Aligned reads to the reference after filtering
* Aligned Filtered Reads %: 100 x Aligned Filtered Reads / Filtered Remaining Reads
* Duplicate Reads #: Number of Duplicates reads (aligned reads having the same 5' alignment positions for both mates in the case of paired-end reads) after filtering
* Duplicate %: 100 x Duplicate Reads / Aligned Filtered Reads
* Final Aligned Reads: Number of Aligned reads without duplicates count (Aligned Filtered Reads - Duplicate Reads)
* Mitchondrial Reads #: Number of reads aligned to chrM
* Mitochondrial %: 100 x Mitochondrial Reads / Aligned Filtered Reads 
