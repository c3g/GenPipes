### Sequencing, Alignment and Methylation Metrics per Sample

General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: Sequencing, Alignment and Methylation Statistics per Sample (**partial table**; [download full table](sequenceAlignmentTable.tsv))

$sequence_alignment_table$

* raw_reads: total number of reads obtained from the sequencer
* trimmed Reads: number of remaining reads after the trimming step
* %_survivalrate: Trimmed reads / Raw reads
* aln_reads: number of aligned reads to the reference
* %_total_aln: Mapped reads / trimmed reads
* DuplicatedReads: number of duplicated read entries providing alternative coordinates
* %aligned_duplicate: DuplicatedReads / aln_reads
* DeduplicatedAlignRreads: aln_reads - DuplicatedReads
* %_UsefulAlignRate: (aln_reads - DuplicatedReads) / raw_reads
* lambdaConversion: C->T conversion rate on the lambda phage
* mean_genomecoverage: total number of aligned reads / genome size
* #_CG_1X: total number of bases with a coverage >= 1x / genome size
* #_CG_10X: total number of bases with a coverage >= 10x / genome size
* #_CG_30X: total number of bases with a coverage >= 30x / genome size
