### FPKM Analysis

#### FPKM Values Generation

The [Cufflinks] program [@cufflinks] is used to assemble aligned RNA-Seq reads into transcripts and to estimate their abundance (FPKM). In RNA-Seq experiments, cDNA fragments are sequenced and mapped back to genes and ideally, individual transcripts. Properly normalized, the RNA-Seq fragment counts can be used as a measure of the relative abundance of transcripts, and Cufflinks measures transcript abundances in **Fragments Per Kilobase of exon per Million fragments mapped (FPKM)**, which is analogous to single-read "RPKM". In paired-end RNA-Seq experiments, fragments are sequenced from both ends, providing two reads for each fragment. To estimate isoform-level abundances, one must assign fragments to individual transcripts, which may be difficult because a read may align to multiple isoforms of the same gene. Cufflinks uses a statistical model of paired-end sequencing experiments to derive a likelihood for the abundances of a set of transcripts given a set of fragments. Once transcripts are assembled and their corresponding FPKM estimated, these transcripts are annotated with the known reference set of transcripts obtained from the Ensembl database.

### Metrics and Exploratory Analysis

We use several metrics and exploratory analysis in order to control the data quality and to verify the biological reliability of the data.

#### FPKM Metrics

The pairwise sample correlation analysis controls the general transcripts expression consistency between samples. It can also check sample mix-up or error in name assignment. Thus, samples belonging to the same design group/condition are expected to show higher level of correlation.

Table: Pairwise Pearson's correlation value per sample ([download full table](corrMatrixSpearman.txt))

$corr_matrix_spearman_table$

Saturation plots show if there is enough sequencing depth to saturate gene expression at various ranges of expression. In RNA-Seq experiments, saturation would be reached when an increment in the number of reads does not result in additional true expressed transcripts being detected. The precision of any sample statistics (FPKM or RPKM) is affected by sample size (sequencing depth). The saturation analysis will resample a series of subsets from total RNA reads and then calculate RPKM values using each subset. By doing this, we are able to check if the current sequencing depth was saturated or not (or if the RPKM values were stable or not). For each sample we estimate the Percent Relative Error (PRE). The PRE measures how the RPKM estimated from a subset of reads deviates from real expression levels) and the median RPKM of the set of transcripts. Saturation plots are generated independently for four different set of transcripts: high, intermediate, moderate and low expressed transcripts (corresponding to quartiles Q1 to Q4 of median RPKM).

#### Exploratory Analysis 

The exploratory analysis is based on the use of gene read counts per million (CPM) and transcripts RPKM values. The main goal of this type of approach is to detect the possible presence of outlier samples and to explore the homogeneity of biological replicates. Thus, if no outlier are detected and replicates are fairly homogenous, the results of following gene and transcript analysis will be reinforced.

First, we analyzed how different samples are connected to each other using principal component analysis (PCA) on gene expression data (log2CPM) without prior statistical filtering. The plot shows the projection of the samples onto the two-dimensional space spanned by the first and second principal component. These are the orthogonal directions in which the data exhibits the largest and second-largest variability. These two components are usually sufficient to differentiate groups of samples describing the principal conditions of the analysis design.

![PCA (first two components) of the gene log2CPM values ([download high-resolution image](pca_log2CPM.pdf))](pca_log2CPM.png)

---


