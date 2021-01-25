### IHEC metrics

Those metrics are calculated according to IHEC standards, following the IHEC [documentation](https://github.com/IHEC/ihec-assay-standards/tree/master/qc_metrics/chip-seq/script). The metrics from [spp]\ [@spp] \(ENCODE metrics\) are also calculated following the [documentation](https://github.com/kundajelab/phantompeakqualtools) and included to the IHEC metrics table.

The alignment file per sample is filtered using [SAMtools]\ [@samtools]. All alignments with MAPQ smaller than **5** and all unmapped and duplicates reads are excluded from the resulting file in BAM format.
General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: IHEC metrics per Sample ([download table](IHEC_metrics_AllSamples.tsv))

* genome_assembly: name of genome assembly used for analysis
* ChIP_type: type of peak called
* treat_name: name of treatment sample
* ctl_name: name of Input sample
* treat_raw_reads: total number of reads obtained from the sequencer for treatment sample(s)
* treat_trimmed_reads: number of remaining reads after the trimming step for treatment sample(s)
* treat_trimmed_frac: 100 x treat_trimmed_reads / treat_raw_reads
* treat_mapped_reads: number of aligned reads to the reference for treatment sample(s)
* treat_mapped_frac: 100 x treat_mapped_reads / treat_trimmed_reads
* treat_dup_reads: number of aligned reads having the same 5' alignment positions (for both mates in the case of paired-end reads) for treatment sample(s)
* treat_dup_frac: 100 x treat_duplicated_reads / treat_mapped_reads
* treat_filtered_reads: number of aligned reads to the reference after filtering for treatment sample(s)
* treat_filtered_frac: 100 x treat_final_reads / treat_trimmed_reads
* treat_MT_reads: number of reads aligned to either chromosome chrM or chromosome MT for treatment sample(s)
* treat_Mt_frac: 100 x treat_MTreads / treat_final_reads
* ctl_raw_reads: total number of reads obtained from the sequencer for Input
* ctl_trimmed_reads: number of remaining reads after the trimming step for Input
* ctl_trimmed_frac: 100 x treat_trimmed_reads / treat_raw_reads
* clt_mapped_reads: number of aligned reads to the reference for Input
* ctl_mapped_frac: 100 x treat_mapped_reads / treat_trimmed_reads
* ctl_dup_reads: number of aligned reads having the same 5' alignment positions (for both mates in the case of paired-end reads) for Input
* ctl_dup_frac: 100 x treat_duplicated_reads / treat_mapped_reads
* ctl_filtered_reads: number of aligned reads to the reference after filtering for Input
* ctl_filtered_frac: 100 x treat_final_reads / treat_trimmed_reads
* ctl_MT_reads: number of reads aligned to either chromosome chrM or chromosome MT for Input
* ctl_Mt_frac: 100 x treat_MTreads / treat_final_reads
* nmb_peaks: number of peaks called
* reads_in_peaks: number of reads in peaks
* frip: Fraction of Reads In Peaks reads_under_peaks / treat_final_reads
* treat_nsc: Normalized Strand cross-correlation coefficient for treatment sample(s)
* ctl_nsc: Normalized Strand cross-correlation coefficient for Input
* treat_rsc: Relative Strand Cross correlation Coefficient for treatment sample(s)
* ctl_rsc: Relative Strand Cross correlation Coefficient for Input
* treat_Quality: quality tag based on thresholded Relative Strand Cross correlation Coefficient for treatment sample
* ctl_Quality: quality tag based on thresholded Relative Strand Cross correlation Coefficient for Input
* singletons: number of singletons for paired-end data sets
* js_dist: Jensen-Shannon distance (only if Input provided)
* chance_div: CHip-seq ANalytics and Confidence Estimation (only if Input provided)
