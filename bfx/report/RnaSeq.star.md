### Read Alignment to a Reference Genome

The filtered reads are aligned to a reference genome. The genome used in this analysis is **_$scientific_name$_** assembly **$assembly$**. Each readset is aligned using [STAR]\ [@star] which creates a Binary Alignment Map file (.bam). Then, all readset BAM files from the same sample are merged into a single global BAM file using [Picard].
