### Sequencing, Alignment and Methylation Metrics per Sample

General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: Sequencing, Alignment and Methylation Statistics per Sample (**partial table**; [download full table](IHEC.sampleMetricsTable.tsv))

$sequence_alignment_table$

* raw_reads: total number of reads obtained from the sequencer
* trimmed_reads: number of remaining reads after the trimming step
* %_survival_rate: trimmed_reads / raw_reads * 100
* aligned_reads: number of aligned reads to the reference
* %_mapping_efficiency: aligned_reads / trimmed_reads * 100
* duplicated_reads: number of duplicated read entries providing the same mapping coordinates (due to PCR duplicates)
* %_duplication_rate: duplicated_reads / aligned_reads * 100 
* deduplicated_aligned_reads: aligned_reads - duplicated_reads
* %_useful_aligned_rate: deduplicated_aligned_reads / raw_reads * 100 
* %_proportion_unique_filtered_reads_MAPQ>10: deduplicated_aligned_reads with mapping quality > 10 / trimmed_reads * 100
* GC_bias: the Pearson correlation between coverage values and GC content for 1000 bins of 100 base pair across genome
* %_pUC19_methylation_rate: 100 - C->T conversion rate on pUC19 * 100 
* %_lambda_conversion_rate: C->T conversion rate on the lambda phage * 100 
* %_human_conversion_rate: estimation of non-CpG methylation * 100 
* estimated_average_genome_coverage: aligned_reads / genome size
* median_CpG_coverage: median read coverage for on-target CpGs
* #_CG_1X: total number of CpGs with a coverage >= 1x
* #_CG_10X: total number of CpGs with a coverage >= 10x
* #_CG_30X: total number of CpGs with a coverage >= 30x

