#### Aligned Reads Base Quality Filtering

This step allows the filtering of aligned reads based on the alignment quality of the BAM file. The alignment file per sample is filtered using [SAMtools]\ [@samtools]. All alignments with MAPQ smaller than **$min_mapq$** and samflag 4 (read unmapped) are excluded from the resulting file in BAM format.
