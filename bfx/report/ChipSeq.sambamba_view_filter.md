#### Aligned Reads Base Quality Filtering

This step allows the filtering of aligned reads based on the alignment quality of the BAM file. The alignment file per sample is filtered using [Sambamba]\ [@sambamba]. All alignments with MAPQ smaller than **$min_mapq$** are excluded from the resulting file in BAM format as long as unmapped reads.
