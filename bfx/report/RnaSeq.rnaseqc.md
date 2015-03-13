### Trimming and Alignment Metrics per Sample

General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: Trimming and Alignment Statistics per Sample (**partial table**; [download full table](trimAlignmentTable.tsv))

$trim_alignment_table$

* Raw Reads: total number of reads obtained from the sequencer
* Surviving Reads: number of remaining reads after the trimming step
* %: Surviving Reads / Raw Reads
* Aligned Reads: number of aligned reads to the reference
* %: Aligned reads / Surviving reads
* Alternative Alignments: number of duplicate read entries providing alternative coordinates
* %: Alternative Alignments / Aligned Reads
* rRNA Reads: number of reads aligning to rRNA regions as defined in the transcript model definition
* %: rRNA Reads / Surviving Reads
* Coverage: mean coverage = number of bp aligning to reference / total number of bp in the reference
* Exonic Rate: fraction mapping reads within exons
* Genes: number of Genes with at least 5 reads

[Additional metrics can be found in the original RNAseqQC report available here](reportRNAseqQC.zip)
