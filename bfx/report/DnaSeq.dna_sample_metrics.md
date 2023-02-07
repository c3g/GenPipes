### Sequencing and Alignment Metrics per Sample

General summary statistics are provided per sample. Sample readsets are merged for clarity.

Table: Sequencing and Alignment Statistics per Sample (**partial table**; [download full table](sequenceAlignmentTable.tsv))

$sequence_alignment_table$

* Raw Reads: total number of reads obtained from the sequencer
* Surviving Reads: number of remaining reads after the trimming step
* Surviving %: Surviving reads / Raw reads
* Mapped reads: number of aligned reads to the reference
* Mapped %: Mapped reads / Surviving reads
* Not Duplicate Reads: number of not duplicated read entries
* Duplicate Reads: number of duplicated read entries providing alternative coordinates
* Duplicate %: Duplicate / Mapped reads
* Pair Orientation: library paired-end read design
* Mean Insert Size: mean distance between the left most base position of the read1 and the right most base position of the read 2
* Standard Deviation: standard deviation of distance between the left most base position of the read1 and the right most base position of the read 2
* Mean Coverage: total number of aligned reads / (genome size or capture region size)
* %_bases_above_10: total number of bases with a coverage >= 10x / (genome size or capture region size)
* %_bases_above_25: total number of bases with a coverage >= 25x / (genome size or capture region size)
* %_bases_above_50: total number of bases with a coverage >= 50x / (genome size or capture region size)
* %_bases_above_75: total number of bases with a coverage >= 75x / (genome size or capture region size)
* %_bases_above_100: total number of bases with a coverage >= 100x / (genome size or capture region size)
* %_bases_above_500: total number of bases with a coverage >= 500x / (genome size or capture region size)

