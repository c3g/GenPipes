### Read discarding and Chimera detection

Reads with at least **$sequence_max_n$** 'N' base are discarded then chimera detection is applied to the merged reads using [Uchime]\ [@uchime]. Uchime performs reference based detection. The **$chimera_db$**\ [@$chimera_ref$] database is used for reference based detection.

Table: Chimera Detection Statistics per Readset ([download full table](uchimeReadsetTable.tsv))

$uchime_readset_table$

* Total Merged $read_type$ Reads #: number of Remaining $read_type$ Reads after the merging step
* Total Chimera filtered out $read_type$ Reads #: number of Remaining $read_type$ Reads after the chimera removal step
* Total Chimera filtered out $read_type$ Reads %: percentage of Total Chimera filtered out $read_type$ Reads / Total Merged $read_type$ Reads
