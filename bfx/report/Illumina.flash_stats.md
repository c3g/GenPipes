### Paired end read Merging

Paired end reads are merged with a minimum overlap of **$min_overlap$** and a maximum overlap of **$max_overlap$**. Merging are performed using [FLASh]\ [@flash].

Table: Merging Statistics per Readset ([download full table](mergeReadsetTable.tsv))

$merge_readset_table$

* Trim $read_type$ Reads #: number of $read_type$ Reads obtained after the trimming step
* Merged $read_type$ Reads #: number of Remaining $read_type$ Reads after the merging step
* Merged $read_type$ Reads %: percentage of Merging $read_type$ Reads / Trim $read_type$ Reads
