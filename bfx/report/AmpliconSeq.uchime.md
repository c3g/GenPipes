### Chimera detection

Chimera detection is applied to the merged reads using [Uchime]\ [@uchime]. Uchime performs both de novo (abundance based) chimera and reference based detection A custom database is used for reference based detection.

Table: Chimera detection Statistics per Readset ([download full table](uchimeReadsetTable.tsv))

$uchime_readset_table$

* Total Merged $read_type$ Reads #: number of Remaining $read_type$ Reads after the merging step
* Total Chimera filtered out $read_type$ Reads #: number of Remaining $read_type$ Reads after the chimera removal step
* Total Chimera filtered out $read_type$ Reads %: percentage of Total Chimera filtered out $read_type$ Reads / Total Merged $read_type$ Reads
