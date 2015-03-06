### Read Alignment to a Reference Genome

The filtered reads are aligned to a reference genome. The genome used in this analysis is **_$scientific_name$_** assembly **$assembly$**. Each readset is aligned using [BWA]\ [@bwa] which creates a Binary Alignment Map file (.bam). Then, all readset BAM files from the same sample are merged into a single global BAM file using [Picard].

BWA is a fast light-weighted tool which aligns relatively short sequences (queries) to a sequence database (target), such as the human reference genome. It's based on the Burrows-Wheeler Transform (BWT). BWA is designed for short queries up to ~200 bp with low error rate (< 3%). It does gapped global alignment, supports paired-end reads, and is one of the fastest short read alignment algorithms to date while also visiting suboptimal hits.
