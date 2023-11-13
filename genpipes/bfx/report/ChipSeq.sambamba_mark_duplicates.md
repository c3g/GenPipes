#### Marking Duplicates

Aligned reads are duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score) will be marked as a duplicate in the BAM file. Duplicates reads will be excluded in the subsequent analysis. Marking duplicates is performed using [Sambamba].
