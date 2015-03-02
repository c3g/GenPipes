### Read Trimming and Clipping of Adapters

Reads are trimmed from the 3' end to have a phred score of at least 30. Illumina sequencing adapters are removed from the reads, and all reads are required to have a length of at least 50 bp. Trimming and clipping are performed using [Trimmomatic] ref. [@trimmomatic].

Table: Trimming Statistics by Readset ([download](trimming.stats))

Sample | Readset | Raw Read # | Surviving Read # | Surviving Read %
-----|-----|-----:|-----:|-----:
