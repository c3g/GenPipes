#### Removing Duplicates

Aligned reads are duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score) will be considered as a duplicate in the BAM file. Duplicates reads will be removed from the bam file and excluded in the subsequent analysis. Removing duplicates is performed using [Picard].
