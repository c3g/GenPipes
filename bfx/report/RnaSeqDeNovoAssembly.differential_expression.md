### Differential Expression Analysis - Methods

#### Differential Analysis Design

The experimental design resulted from a discussion with the client. The designs used in differential analysis are presented in the following table, which contains the sample names as well as the sample group membership per design. For each experimental design (column name), three conditions/groups are possible: 0, 1 and 2. If a sample is assigned 0, it is not included in a particular analysis. If a sample is assigned 1, the sample is considered as a member of the control group. If a sample is assigned 2, the sample is considered as a member of the test/case group.

Table: Sample Names and Experimental Designs (**partial table**; [download full table](design.tsv))

$design_table$


#### Differential Gene Analysis Description

A primary task in the analysis of RNA-Seq data is the detection of differentially expressed genes. For this purpose, count data from RNA-Seq should be obtained for each non-overlapping gene. Read counts are found to be (to good approximation) linearly related to the abundance of the target transcript [@quantifying_rnaseq]. If reads were independently sampled from a population with given, fixed fractions of genes, the read counts would follow a multinomial distribution, which can be approximated by the Poisson distribution. Thus, we can use statistical testing to decide whether, for a given gene, an observed difference in read counts is significant, that is, whether it is greater than what would be expected just due to natural random variation. The differential gene expression analysis is done using [DESeq]\ [@deseq] and [edgeR]\ [@edger] R Bioconductor packages

Gene abundance estimation is performed using [RSEM]\ [@rsem] and is represented as a table which reports, for each sample (columns), the number of reads mapped to a given gene (rows).

Table: Matrix of Raw Read Counts per Gene per Sample (**partial table**; [download full table](rawCountMatrix.csv))

$raw_count_matrix_table$


### Differential Expression Analysis - Results

The following sections provide the results of the differential gene expresssion per design.

* Gene: component id
* Symbol: BLAST protein database id if available, component id otherwise
* log_FC: log2 Fold Change of gene level expression
* log_CPM: log2 Counts Per Million of gene level expression
* deseq.p-value: DESeq nominal p-value
* deseq.adj.pvalue: DESeq False Discovery Rate (FDR) adjusted p-value
* edger.p-value: edgeR nominal p-value
* edger.adj.pvalue: edgeR False Discovery Rate (FDR) adjusted p-value

