### IHEC metrics

Those metrics are calculated according to IHEC standards, following the IHEC [documentation](https://github.com/IHEC/ihec-assay-standards/tree/master/qc_metrics/chip-seq/script). The metrics from [spp]\ [@spp] \(ENCODE metrics\) are also calculated following the [documentation](https://github.com/kundajelab/phantompeakqualtools) and included to the IHEC metrics table.

The alignment file per sample is filtered using [SAMtools]\ [@samtools]. All alignments with MAPQ smaller than **5** and all unmapped and duplicates reads are excluded from the resulting file in BAM format.
General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: IHEC metrics per Sample ([download table](IHEC_chipseq_metrics_AllSamples.tsv))

* Sample_Name: Name of sample
* Mark_Name: Name of Histone mark
* ChIP_type: Type of peak called
* Genome_Assembly: Genome Assembly used for analysis
* Raw_Reads: Total Number of reads obtained from the sequencer
* Trimmed_Reads: Number of remaining reads after Trimming
* Trimmed_Reads_Fraction: 100 x Trimmed_Reads / Raw_Reads
* Mapped_Reads: Number of Aligned reads to the reference after Trimming
* Mapped_Reads_Fraction: 100 x Mapped_Reads / Trimmed_Reads
* Duplicates_Reads: Number of Duplicates reads ( aligned reads to the reference having the same 5' alignment positions for both mates in the case of paired-end reads)
* Duplicates_Reads_Fraction: 100 x Duplicates_Reads / Mapped_Reads
* Filtered_Reads: Number of Aligned reads to the reference after filtering
* Filtered_Reads_Fraction: 100 x Filtered_Reads / Trimmed_Reads
* Mitochondrial_Reads: Number of reads Aligned to either chromosome chrM or chromosome MT
* Mitochondrial_Reads_Fraction: 100 x Mitochondrial_Reads / Filtered_Reads
* Singletons_Reads: Number of Singletons for Paired-End data sets
* Nbr_Peaks: Number of peaks called
* Reads_in_Peaks: Number of reads in peaks
* frip: Fraction of Reads In Peaks Reads_in_Peaks / Filtered_Reads
* Mark_NSC: Normalized Strand Cross-Correlation coefficient
* Mark_RSC: Relative Strand Cross-Correlation coefficient
* Mark_Quality: Quality Tag based on thresholded RSC coefficient
* JS_Distance: Jensen-Shannon distance (only if Input provided)
* Chance_Divergence: ChIP-Seq Analytics and Confidence Estimation (only if Input provided)