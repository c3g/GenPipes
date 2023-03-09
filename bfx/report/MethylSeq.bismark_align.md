### Read Alignment to a Reference Genome

The filtered reads are aligned to a reference genome. The genome used in this analysis is **_$scientific_name$_** assembly **$assembly$**. Each readset is aligned using [Bismark]\ [@bismark] which creates a Binary Alignment Map file (.bam). Then all readgroup are added to all readset BAM files using [Picard]. Finally , all readset BAM files from the same sample are merged into a single global BAM file, also using [Picard].
