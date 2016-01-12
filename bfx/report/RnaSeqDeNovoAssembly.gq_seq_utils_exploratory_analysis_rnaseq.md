#### Exploratory Analysis 

The exploratory analysis is based on the use  of [RSEM]\ [@rsem] gene read counts per million (CPM). The main goal of this type of approach is to detect the possible presence of outlier samples and to explore the homogeneity of biological replicates. Thus, if no outlier are detected and replicates are fairly homogenous, the results of following gene and transcript analysis will be reinforced.

First, we analyze how different samples are connected to each other using principal component analysis (PCA) on gene expression data (log2CPM) without prior statistical filtering. The plot shows the projection of the samples onto the two-dimensional space spanned by the first and second principal component. These are the orthogonal directions in which the data exhibits the largest and second-largest variability. These two components are usually sufficient to differentiate groups of samples describing the principal conditions of the analysis design.

Secondly, we use RPKM values of the entire set of transcripts to estimate the correlation distance between each sample. A hierarchical clustering of these distances (Ward approach) is realized to visualize the transcript expression divergences between samples.

Finally, genes with the most variable expression data (log2CPM standard deviation) are used to vizualize how the most variable genes are expressed among samples. The heatmap plot indicates if it exists specific patterns of variation among these genes that could be used to discriminate between groups of samples.

