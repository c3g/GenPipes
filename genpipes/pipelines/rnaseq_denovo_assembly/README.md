<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [RNA-Seq De Novo Assembly Pipeline](#rna-seq-de-novo-assembly-pipeline)
  - [Usage](#usage)
  - [picard_sam_to_fastq](#picard_sam_to_fastq)
  - [trimmomatic](#trimmomatic)
  - [merge_trimmomatic_stats](#merge_trimmomatic_stats)
  - [insilico_read_normalization_readsets](#insilico_read_normalization_readsets)
  - [insilico_read_normalization_all](#insilico_read_normalization_all)
  - [trinity](#trinity)
  - [exonerate_fastasplit](#exonerate_fastasplit)
  - [blastx_trinity_uniprot](#blastx_trinity_uniprot)
  - [blastx_trinity_uniprot_merge](#blastx_trinity_uniprot_merge)
  - [transdecoder](#transdecoder)
  - [hmmer](#hmmer)
  - [infernal_transcriptome](#infernal_transcriptome)
  - [blastp_transdecoder_uniprot](#blastp_transdecoder_uniprot)
  - [signalp](#signalp)
  - [tmhmm](#tmhmm)
  - [trinotate](#trinotate)
  - [align_and_estimate_abundance_prep_reference](#align_and_estimate_abundance_prep_reference)
  - [align_and_estimate_abundance](#align_and_estimate_abundance)
  - [gq_seq_utils_exploratory_analysis_rnaseq_denovo](#gq_seq_utils_exploratory_analysis_rnaseq_denovo)
  - [differential_expression](#differential_expression)
  - [filter_annotated_components](#filter_annotated_components)
  - [gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered](#gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered)
  - [differential_expression_filtered](#differential_expression_filtered)
  - [multiqc](#multiqc)
  - [merge_fastq](#merge_fastq)
  - [seq2fun](#seq2fun)
  - [differential_expression_seq2fun](#differential_expression_seq2fun)
  - [pathway_enrichment_seq2fun](#pathway_enrichment_seq2fun)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

RNA-Seq De Novo Assembly Pipeline
=================================

The standard Genpipes RNA-Seq De Novo Assembly pipeline now has two protocols:
The [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
software suite to reconstruct transcriptomes from RNA-Seq data or the Seq2Fun software suite to generate several
informative outputs including gene abundance tables, pathway and species hit tables which
are required for the [Networkanalyst] (https://www.networkanalyst.ca/home.xhtml). Only RNA-Seq data is used for
both of the protocols and any reference genome or transcriptome is not required.

The trinity De Novo Assembly pipeline, selected using the "-t trinity" parameter starts by trimming reads
with [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
Then reads are normalized in order to reduce memory requirement and decrease assembly runtime, using the Trinity
normalization utility inspired by the [Diginorm](http://arxiv.org/abs/1203.4802) algorithm.

Then, the transcriptome is assembled on normalized reads using the Trinity assembler. Trinity creates
a Trinity.fasta file with a list of contigs representing the transcriptome isoforms. Those transcripts
are grouped in components mostly representing genes.

Components and transcripts are functionally annotated using the [Trinotate](https://github.com/Trinotate/Trinotate/wiki) suite.

Gene abundance estimation for each sample has been performed using [RSEM](https://deweylab.github.io/RSEM/)
(RNA-Seq by Expectation-Maximization). Differential gene expression analysis is performed using
[DESeq](http://genomebiology.com/2010/11/10/R106) and [edgeR](http://bioinformatics.oxfordjournals.org/content/26/1/139/) R Bioconductor packages.

The DESeq and edgeR methods model **count data** by a negative binomial distribution. The parameters of
the distribution (mean and dispersion) are estimated from the data, i.e. from the read counts in the input files.
Both methods compute a measure of read abundance, i.e. expression level (called *base mean* or
*mean of normalized counts* in DESeq, and *concentration* in edgeR) for each gene and apply a hypothesis test
to each gene to evaluate differential expression. In particular, both methods determine a p-value and
a log2 fold change (in expression level) for each gene. The Log2 FC of EdgeR is reported in the differential gene
results file, one file per design.

The log2fold change is the logarithm (to basis 2) of the fold change condition from condition A to B
(mutation or treatment are the most common conditions). A "fold change" between conditions A and B at a gene
or transcript is normally computed as the ratio at gene or transcript of the base mean of scaled counts
for condition B to the base mean of scaled counts for condition A. Counts are scaled by a size factor in
a step called normalization (if the counts of non-differentially expressed genes in one sample are, on average,
twice as high as in another,  the size factor for the first sample should be twice that of the other sample).
Each column of the count table is then divided by the size factor for this column and the count values
are brought to a common scale, making them comparable. See the [EdgeR vignette](http://www.bioconductor.org/packages/2.12/bioc/vignettes/edgeR/inst/doc/edgeR.pdf) for additional information on normalization approaches used in the pipeline.

The differential gene analysis is followed by a Gene Ontology (GO) enrichment analysis.
This analysis use the [goseq approach](http://bioconductor.org/packages/release/bioc/html/goseq.html).
The goseq is based on the use of non-native GO terms resulting from trinotate annotations (see details in the section 5 of
[the corresponding vignette](http://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf).

Thus a high quality contigs assembly is created by extracting all transcripts having a functionnal annotation as defined by trinotate,
the Top BLASTX hit and TmHMM annotations are used by default.

Finally, different exploratory data analysis (EDA) techniques are applied to filtered isoforms expression levels.
Main goals of expression level EDA are the detection of outliers, potential mislabeling,  to explore the homogeneity
of biological replicates and  to appreciate the global effects of the different experimental variables.

An HTML summary report is automatically generated by the pipeline. Various Quality Control (QC) 
summary statistics are included in the report and can be explored interactively.
    
The Seq2Fun De Novo Assembly pipeline, selected using the "-t seq2fun" parameter directly starts with Seq2Fun
software suit from fastq files.
    
Usage
-----

```
#!text
usage: genpipes rnaseq_denovo_assembly [-h] [--clean] -c CONFIG [CONFIG ...]
                                       [--container {wrapper, singularity} <IMAGE PATH>]
                                       [-f]
                                       [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                                       [--genpipes_file GENPIPES_FILE]
                                       [-j {pbs,batch,daemon,slurm}]
                                       [--json-pt]
                                       [-l {debug,info,warning,error,critical}]
                                       [-o OUTPUT_DIR] [--sanity-check]
                                       [-s STEPS] [--wrap [WRAP]]
                                       -r READSETS_FILE [-d DESIGN_FILE] [-v]
                                       [-t {trinity,seq2fun}] [-b BATCH]

For more documentation, visit our website: https://genpipes.readthedocs.io

options:
  -h, --help            show this help message and exit
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -c, --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity
                        image path
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --force_mem_per_cpu FORCE_MEM_PER_CPU
                        Take the mem input in the ini file and force to have a
                        minimum of mem_per_cpu by correcting the number of cpu
                        (default: None)
  --genpipes_file, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -j, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
  -l, --log {debug,info,warning,error,critical}
                        log level (default: info)
  -o, --output-dir OUTPUT_DIR
                        output directory (default: current)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  -s, --steps STEPS     step range e.g. '1-5', '3,6,7', '2,4-8'
  --wrap [WRAP]         Path to the genpipes cvmfs wrapper script. Default is 
                        genpipes/ressources/container/bin/container_wrapper.sh
                        . This is a convenience option for using genpipes in a
                        container
  -r, --readsets READSETS_FILE
                        readset file
  -d, --design DESIGN_FILE
                        design file
  -v, --version         show the version information and exit
  -t, --type {trinity,seq2fun}
                        RNAseq analysis type
  -b, --batch BATCH     batch file (to peform batch effect correction

Steps:

Protocol trinity
1 picard_sam_to_fastq
2 trimmomatic
3 merge_trimmomatic_stats
4 insilico_read_normalization_readsets
5 insilico_read_normalization_all
6 trinity
7 exonerate_fastasplit
8 blastx_trinity_uniprot
9 blastx_trinity_uniprot_merge
10 transdecoder
11 hmmer
12 infernal_transcriptome
13 blastp_transdecoder_uniprot
14 signalp
15 tmhmm
16 trinotate
17 align_and_estimate_abundance_prep_reference
18 align_and_estimate_abundance
19 gq_seq_utils_exploratory_analysis_rnaseq_denovo
20 differential_expression
21 filter_annotated_components
22 gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered
23 differential_expression_filtered
24 multiqc

Protocol seq2fun
1 picard_sam_to_fastq
2 merge_fastq
3 seq2fun
4 differential_expression_seq2fun
5 pathway_enrichment_seq2fun
```

picard_sam_to_fastq 
-------------------
 
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

trimmomatic 
-----------
 
Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', parameter 'adapter_fasta'),
it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
only Adapter1 is used and left unchanged.

This step takes as input files:
1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

merge_trimmomatic_stats 
-----------------------
 
The trim statistics per readset are merged at this step.

insilico_read_normalization_readsets 
------------------------------------
 
Normalize each readset, using the Trinity normalization utility.

insilico_read_normalization_all 
-------------------------------
 
Merge all normalized readsets together and normalize the result, using the Trinity normalization utility.

trinity 
-------
 
Create a de novo assembly from normalized readsets using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki).

exonerate_fastasplit 
--------------------
 
Split the Trinity assembly FASTA into chunks for further parallel BLAST annotations.

blastx_trinity_uniprot 
----------------------
 
Annotate Trinity FASTA chunks with Swiss-Prot (default), UniRef or a custom database using [blastx](http://blast.ncbi.nlm.nih.gov/).

blastx_trinity_uniprot_merge 
----------------------------
 
Merge blastx Swiss-Prot chunks results.

transdecoder 
------------
 
Identifies candidate coding regions within transcript sequences using [Transdecoder](http://transdecoder.github.io/).

hmmer 
-----
 
Identifies protein domains using [HMMR](http://hmmer.janelia.org/).

infernal_transcriptome 
----------------------
 
Identify structural RNAs using cmscan function from [infernal](http://eddylab.org/infernal)
Run in parallel, using chunks created by exonerate. 

blastp_transdecoder_uniprot 
---------------------------
 
Search Transdecoder-predicted coding regions for sequence homologies on UniProt using [blastp](http://blast.ncbi.nlm.nih.gov/).

signalp 
-------
 
Predict signal peptides using [SignalP](https://services.healthtech.dtu.dk/services/SignalP-6.0/).

tmhmm 
-----
 
Predict transmembrane regions using [TMHMM](https://services.healthtech.dtu.dk/services/TMHMM-2.0/).

trinotate 
---------
 
Perform transcriptome functional annotation and analysis using [Trinotate](https://github.com/Trinotate/Trinotate/wiki).
All functional annotation data is integrated into a SQLite database and a whole annotation report is created.

align_and_estimate_abundance_prep_reference 
-------------------------------------------
 
Index Trinity FASTA file for further abundance estimation using [Trinity align_and_estimate_abundance.pl utility](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification).

align_and_estimate_abundance 
----------------------------
 
Estimate transcript abundance using [RSEM](http://deweylab.github.io/RSEM/) via
[Trinity align_and_estimate_abundance.pl utility](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification).

gq_seq_utils_exploratory_analysis_rnaseq_denovo 
-----------------------------------------------
 
Exploratory analysis using the gqSeqUtils R package.

differential_expression 
-----------------------
 
Performs differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
Merge the results of the analysis in a single csv file. Also, performs Gene Ontology analysis for RNA-Seq denovo Assembly using the Bioconductor's R package [goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html).
Generates GO annotations for differential genes and isoforms expression analysis, based on associated GOTERMS generated by trinotate.

filter_annotated_components 
---------------------------
 
Filter high quality contigs based on values in trinotate annotations. Recreate a high quality contigs fasta file and run Assembly statistics using the gqSeqUtils R package.

gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered 
--------------------------------------------------------
 
Exploratory analysis using the gqSeqUtils R package using a subset of filtered transcripts.

differential_expression_filtered 
--------------------------------
 
Differential Expression and GOSEQ analysis based on filtered transcripts and genes.

multiqc 
-------
 
Aggregate results from bioinformatics analyses across many samples into a single report.
MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).

merge_fastq 
-----------
 
This step is performed to merge fastq files if multiple readset files for one sample is present

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

seq2fun 
-------
 
seq2fun

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

This step perform seq2fun analysis and generates output files including KO abundance table and KO mapped fastq files
(https://www.seq2fun.ca/manual.xhtml#sect4) and (https://www.seq2fun.ca/manual.xhtml#sect20)

For each contrast different folders and all the files for that particular contrast are
generated. Therefore, only pairwise comparisons are possible
(treatment and controls will be added according to the 1 and 2 in the design file).

differential_expression_seq2fun 
-------------------------------
 
Performs differential gene expression analysis using [DESEQ2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
Merge the results of the analysis in a single csv file.

pathway_enrichment_seq2fun 
--------------------------
 

seq2fun pathway analysis using fgsea (https://bioconductor.org/packages/release/bioc/html/fgsea.html)
 and user provide universal pathway list as KEGG map ID. The differential KO expression results obtained
 from edgeR will be using as the input for the pathway enrichment analysis

