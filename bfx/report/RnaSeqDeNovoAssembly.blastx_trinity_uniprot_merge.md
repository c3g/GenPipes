### BLAST Annotation

Each transcript has been aligned against the **_$blast_db$_** protein database using `blastx` program from the NCBI [BLAST] family.

For each transcript, the best BLAST hit is used to annotate its associated component/gene in further Differential Expression Results.

[The full BLAST results file is available here](blastx_Trinity_$blast_db$.tsv.zip)

Table: BLAST Tabular Output Default Columns Description 

|Column_number|Column_id|Column_description|
|:------------ |---------------:| -----:|
|1.|qseqid|query|(e.g.,gene) sequence id |
|2.|sseqid|subject|(e.g., reference genome) sequence id|
|3.|pident|percentage of identical matches|
|4.|length|alignment length|
|5.|mismatch|number of mismatches|
|6.|gapopen|number of gap|openings|
|7.|qstart|start of alignment in query|
|8.|qend|end of alignment in query|
|9.|sstart|start of alignment in subject|
|10.|send|end of alignment in subject|
|11.|evalue|expect value|
|12.|bitscore|bit score|
