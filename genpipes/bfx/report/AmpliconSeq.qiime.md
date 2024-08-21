### OTU generation

Taxomomic affiliation of the **$amplicon_type$** data has been studied using [Qiime]\ [@qiime].

#### OTU picking

Demultiplexed and quality filtered sequences from pre-processing step are clustered into OTUs using [VSEARCH]\ [@vsearch]. An OTU (Operational Taxonomic Unit) is formed based on sequence identity: the identity threshold defined is **$similarity$**.

After the clustering step, a representative sequence is picked for each OTU. Each OTU is therefore represented by a single sequence.

#### Classification

A taxonomic identity is assigned to each representative sequence. The database used for **$amplicon_type$** data is **$amplicon_db$**. [Uclust]\ [@uclust] has been used for taxonomic assignation.

#### Alignment and Phylogenetic tree

Multiple alignment of the representative OTU sequences is generated with [PyNAST]\ [@pynast]. It aligns the sequences to a template alignment of reference  **$amplicon_type$** sequences. Then a phylogenetic tree is built in order to study the relationship between the sequences and the samples by computing UniFrac distances. Phylogenetic tree is generated with [FastTree]\ [@fasttree].

