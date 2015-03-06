### Read Trimming and Clipping of Adapters

Reads are trimmed from the 3' end to have a phred score of at least 30. Illumina sequencing adapters are removed from the reads, and all reads are required to have a length of at least 50 bp. Trimming and clipping are performed using [Trimmomatic]\ [@trimmomatic].

Table: Trimming Statistics per Readset ([download full table](trimReadsetTable.tsv))

$trim_readset_table$

* Raw $read_type$ Reads #: number of $read_type$ Reads obtained from the sequencer
* Surviving $read_type$ Reads #: number of Remaining $read_type$ Reads after the trimming step
* Surviving $read_type$ Reads %: percentage of Surviving $read_type$ Reads / Raw $read_type$ Reads
