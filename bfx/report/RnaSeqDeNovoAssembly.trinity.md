### De Novo Assembly

The transcriptome has been assembled on normalized reads using the [Trinity] assembler [@trinity]. Briefly, the assembly process consists of 3 steps:

* **_Inchworm_** assembles the RNA-Seq data into the unique sequences of transcripts, often generating full-length transcripts for a dominant isoform, but then reports just the unique portions of alternatively spliced transcripts.
* **_Chrysalis_** clusters the Inchworm contigs into clusters and constructs complete de Bruijn graphs for each cluster. Each cluster represents the full transcriptonal complexity for a given gene (or sets of genes that share sequences in common). Chrysalis then partitions the full read set among these disjoint graphs.
* **_Butterfly_** then processes the individual graphs in parallel, tracing the paths that reads and pairs of reads take within the graph, ultimately reporting full-length transcripts for alternatively spliced isoforms, and teasing apart transcripts that corresponds to paralogous genes.

Trinity has created a `Trinity.fasta` file with a list of contigs representing the transcriptome isoforms. Those transcripts are grouped in components loosely representing genes.
Transcript names are prefixed by the component/gene name e.g. transcripts `c115_g5_i1` and `c115_g5_i2` are derived from the same isolated de Bruijn graph and therefore share the same component/gene number `c115_g5`.

[The full Trinity assembly FASTA file is available here](Trinity.fasta.zip)


Table: Trinity De Novo Assembly Metrics (**partial table**; [download full table](Trinity.stats.pdf))

Description|Value
-----|----:
$assembly_table$

![Transcript Length Distribution](Trinity.stats.jpg)
