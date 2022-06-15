[TOC]


Amplicon-Seq Pipeline
================



Usage
-----
```
#!text

usage: ampliconseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                      [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                      [--no-json] [--report] [--clean]
                      [-l {debug,info,warning,error,critical}]
                      [--sanity-check]
                      [--container {wrapper, singularity} <IMAGE PATH>]
                      [--genpipes_file GENPIPES_FILE] [-t {qiime,dada2}]
                      [-d DESIGN] [-r READSETS] [-v]

Version: 4.3.0

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  -j {pbs,batch,daemon,slurm}, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
  --report              create 'pandoc' command to merge all job markdown
                        report files in the given step range into HTML, if
                        they exist; if --report is set, --job-scheduler,
                        --force, --clean options and job up-to-date status are
                        ignored (default: false)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -t {qiime,dada2}, --type {qiime,dada2}
                        AmpliconSeq analysis type
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
```
![ampliconseq qiime workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_ampliconseq_qiime.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_ampliconseq_qiime.png)
```
qiime:
1- trimmomatic16S
2- merge_trimmomatic_stats16S
3- flash_pass1
4- ampliconLengthParser
5- flash_pass2
6- merge_flash_stats
7- catenate
8- uchime
9- merge_uchime_stats
10- otu_picking
11- otu_rep_picking
12- otu_assigning
13- otu_table
14- otu_alignment
15- filter_alignment
16- phylogeny
17- qiime_report
18- multiple_rarefaction
19- alpha_diversity
20- collate_alpha
21- sample_rarefaction_plot
22- qiime_report2
23- single_rarefaction
24- css_normalization
25- rarefaction_plot
26- summarize_taxa
27- plot_taxa
28- plot_heatmap
29- krona
30- plot_to_alpha
31- beta_diversity
32- pcoa
33- pcoa_plot
34- plot_to_beta
----
```
![ampliconseq dada2 workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_ampliconseq_dada2.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_ampliconseq_dada2.png)
```
dada2:
1- trimmomatic16S
2- merge_trimmomatic_stats16S
3- flash_pass1
4- ampliconLengthParser
5- flash_pass2
6- merge_flash_stats
7- asva

```

trimmomatic16S
--------------
MiSeq raw reads adapter & primers trimming and basic QC is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
it is used first. Else, Adapter1, Adapter2, Primer1 and Primer2 columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. Sequences are reversed-complemented and swapped.

This step takes as input files:
1. MiSeq paired-End FASTQ files from the readset file

merge_trimmomatic_stats16S
--------------------------
The trim statistics per readset are merged at this step.

flash_pass1
-----------
ampliconLengthParser
--------------------
look at FLASH output to set amplicon lengths input for dada2. As minimum elligible length, a given length needs to have at least 1% of the total number of amplicons

flash_pass2
-----------
merge_flash_stats
-----------------
The paired end merge statistics per readset are merged at this step.

catenate
--------
Catenate all the reads in one file for further analysis.

This step takes as input files:

1. Merged FASTQ files from previous flash step.

uchime
------
Reference based chimera detection is performed using [vsearch](https://github.com/torognes/vsearch)

This step takes as input files:

1. Catenated FASTA file from previous step catenate.

merge_uchime_stats
------------------
The chimeric sequences filtered out statistics per readset are merged at this step.

otu_picking
-----------
The OTU picking step (de novo & close_ref) assigns similar sequences to operational taxonomic units (OTUs) by clustering sequences based on a user-defined similarity threshold. Method per default uses [VSEARCH] (https://github.com/torognes/vsearch) and [Qiime] (http://qiime.org).

This step takes as input file:

1. Catenated and filtered FASTA file from previous step.

otu_rep_picking
---------------
After picking OTUs, this step pick a representative sequence for each OTU.

This step takes as input files:

1. OTU file from previous step
2. Catenated and filtered FASTA file from filter_chimeras step.

otu_assigning
-------------
Given a set of OTUS, this step attempts to assign the taxonomy of each OTU using [Uclust] (http://drive5.com/usearch/manual/uclust_algo.html).

This step takes as input files:

1. OTU representative sequence file from previous step.

otu_table
---------
This step make a consensus OTU table in biom format. It tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU.

This step takes as input files:

1. OTU picking file.
2. Taxonomy assignment for each OTU from the previous step.

otu_alignment
-------------
Align OTU representative sequences using [PyNAST] (http://biocore.github.io/pynast/).

This step takes as input file:

1. OTU representative sequence file.

filter_alignment
----------------
Filter the alignment by removing positions which are gaps in every sequence.

This step takes as input file:

1. Alignment sequence file.

phylogeny
---------
Build a phylogenetic tree from a multiple sequence alignment using [FastTree] (http://www.microbesonline.org/fasttree/).

This step takes as input file:

1. Filtered alignment sequence file from previous step.

qiime_report
------------
1st part report for taxonomic affiliation.

multiple_rarefaction
--------------------
1st step (/4) for rarefaction plot.
Rarefies OTU table by random sampling (without replacement) at different depth in order to perform rarefaction analysis.
You need to provide the minimum/maximum number of sequences per samples and the size of each steps between the min/max of seqs/sample.

This step takes as input files:

1. OTU non rarefied table in biom format.

alpha_diversity
---------------
2nd step (/4) for rarefaction plot.
Calculate alpha diversity on each sample using a variety of alpha diversity metrics (chao1, shannon, observed otus).

This step takes as input files:

1. Multiple OTU rarefied table in biom format from previous step.

collate_alpha
-------------
3rd step (/4) for rarefaction plot.
Merge all the alpha diversity computed in the previous step.

sample_rarefaction_plot
-----------------------
Last step for rarefaction plot.
Plot the rarefaction curve for each sample

qiime_report2
-------------
2nd part report for taxonomic affiliation. Plot rarefaction curve for each sample.

single_rarefaction
------------------
This step is recommended. It subsamples (rarefy) all the samples to an equal number of sequences for further comparaison.
You have to provide the number of sequences to subsample per sample in the configuration file (single_rarefaction_depth).

This step takes as input files:

1. OTU table in biom format.

css_normalization
-----------------
This step is recommended. Alternative method for normalization to rarefaction.
Performs the CSS Matrix normalization.

This step takes as input files:

1. OTU table in biom format.

rarefaction_plot
----------------
Last step for rarefaction plot.
Rarefaction curve for each sample on the same plot.

summarize_taxa
--------------
1st step (/3) for taxonomic affiliation plot.
Summarize information of taxonomic groups within each sample at different taxonomic level.

This step takes as input files:

1. OTU rarefied table in biom format if available.
2. Else, OTU non rarefied table in biom format.

plot_taxa
---------
2nd step (/3) for taxonomic affiliation plot.
Make taxonomy summary bar plots based on taxonomy assignment.

This step takes as input files:

1. Summarized information from previous step.

plot_heatmap
------------
Last step for taxonomic affiliation plot.
Make heatmap at phylum level.

This step takes as input files:

1. Summarized information from previous step.

krona
-----
Plot Krona chart for taxonomic affiliation

plot_to_alpha
-------------
Final report 1st part for the Amplicon-Seq pipeline. Display results (taxonomy, heatmap and alpha diversity).

beta_diversity
--------------
1st step (/3) for 2D PCoA plot.
Calculate beta diversity (pairwise sample dissimilarity) on OTU table. The OTU table has to be normalized.
Only works with >= 4 samples

This step takes as input files:

1. OTU rarefied table in biom format.
2. Tree file.

pcoa
----
2nd step (/3) for 2D PCoA plot.
Compute coordinates pour PCoA

This step takes as input file:

1. Matrix produced in the previous step.

pcoa_plot
---------
Last step for 2D PCoA plot.

This step takes as input file:

1. PCoA from the previous step.

plot_to_beta
------------
Final report's 2nd part for the Amplicon-Seq pipeline. Display results (beta diversity PCoA plots).

asva
----
check for design file (required for PCA plots)


