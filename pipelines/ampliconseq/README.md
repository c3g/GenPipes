[TOC]


Amplicon-Seq Pipeline
================

T


Usage
-----
```
#!text

usage: ampliconseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                      [-o OUTPUT_DIR] [-j {pbs,batch}] [-f] [--report]
                      [--clean] [-l {debug,info,warning,error,critical}]
                      [-r READSETS] [-v]

Version: 2.2.0-beta

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

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
  -j {pbs,batch}, --job-scheduler {pbs,batch}
                        job scheduler type (default: pbs)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
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
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- flash
5- merge_flash_stats
6- catenate
7- uchime
8- merge_uchime_stats
9- otu_ref_picking
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

```
1- picard_sam_to_fastq
----------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

2- trimmomatic
--------------
Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
only Adapter1 is used and left unchanged.

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

3- merge_trimmomatic_stats
--------------------------
The trim statistics per readset are merged at this step.

4- flash
--------
Merge paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/).

5- merge_flash_stats
--------------------
The paired end merge statistics per readset are merged at this step.

6- catenate
-----------
Catenate all the reads in one file for further analysis.

This step takes as input files:

1. Merged FASTQ files from previous step flash. 

7- uchime
---------
Reference based chimera detection is performed using [vsearch](https://github.com/torognes/vsearch)

This step takes as input files:

1. Catenated FASTA file from previous step catenate. 

8- merge_uchime_stats
---------------------
The chimeric sequences filtered out statistics per readset are merged at this step.

9- otu_ref_picking
------------------
The OTU picking step (closed_ref) assigns similar sequences to operational taxonomic units (OTUs) by clustering sequences based on a user-defined similarity threshold. Method per default uses [VSEARCH] (https://github.com/torognes/vsearch) and [Qiime] (http://qiime.org).

This step takes as input file:

1. Catenated and filtered FASTA file from previous step.


10- otu_picking
---------------
The OTU picking step (de novo) assigns similar sequences to operational taxonomic units (OTUs) by clustering sequences based on a user-defined similarity threshold. Method per default uses [VSEARCH] (https://github.com/torognes/vsearch) and [Qiime] (http://qiime.org).

This step takes as input file:

1. Catenated and filtered FASTA file from previous step.


11- otu_rep_picking
-------------------
After picking OTUs, this step pick a representative sequence for each OTU.

This step takes as input files:

1. OTU file from previous step 
2. Catenated and filtered FASTA file from filter_chimeras step.


12- otu_assigning
-----------------
Given a set of OTUS, this step attempts to assign the taxonomy of each OTU using [Uclust] (http://drive5.com/usearch/manual/uclust_algo.html).

This step takes as input files:

1. OTU representative sequence file from previous step.


13- otu_table
-------------
This step make a consensus OTU table in biom format. It tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU. 

This step takes as input files:

1. OTU picking file.
2. Taxonomy assignment for each OTU from the previous step.


14- otu_alignment
-----------------
Align OTU representative sequences using [PyNAST] (http://biocore.github.io/pynast/).

This step takes as input file:

1. OTU representative sequence file.


15- filter_alignment
--------------------
Filter the alignment by removing positions which are gaps in every sequence.

This step takes as input file:

1. Alignment sequence file.


16- phylogeny
-------------
Build a phylogenetic tree from a multiple sequence alignment using [FastTree] (http://www.microbesonline.org/fasttree/).

This step takes as input file:

1. Filtered alignment sequence file from previous step.


17- qiime_report
----------------
1st part report for taxonomic affiliation. 

18- multiple_rarefaction
------------------------
1st step (/4) for rarefaction plot.
Rarefies OTU table by random sampling (without replacement) at different depth in order to perform rarefaction analysis. 
You need to provide the minimum/maximum number of sequences per samples and the size of each steps between the min/max of seqs/sample. 

This step takes as input files:

1. OTU non rarefied table in biom format.


19- alpha_diversity
-------------------
2nd step (/4) for rarefaction plot.
Calculate alpha diversity on each sample using a variety of alpha diversity metrics (chao1, shannon, observed otus). 

This step takes as input files:

1. Multiple OTU rarefied table in biom format from previous step.


20- collate_alpha
-----------------
3rd step (/4) for rarefaction plot.
Merge all the alpha diversity computed in the previous step. 

21- sample_rarefaction_plot
---------------------------
Last step for rarefaction plot.
Plot the rarefaction curve for each sample

22- qiime_report2
-----------------
2nd part report for taxonomic affiliation. Plot rarefaction curve for each sample.

23- single_rarefaction
----------------------
This step is recommended. It subsamples (rarefy) all the samples to an equal number of sequences for further comparaison.
You have to provide the number of sequences to subsample per sample in the configuration file (single_rarefaction_depth).

This step takes as input files:

1. OTU table in biom format.


24- css_normalization
---------------------
This step is recommended. Alternative method for normalization to rarefaction. 
Performs the CSS Matrix normalization.

This step takes as input files:

1. OTU table in biom format.


25- rarefaction_plot
--------------------
Last step for rarefaction plot.
Rarefaction curve for each sample on the same plot. 

26- summarize_taxa
------------------
1st step (/3) for taxonomic affiliation plot.
Summarize information of taxonomic groups within each sample at different taxonomic level. 

This step takes as input files:

1. OTU rarefied table in biom format if available.
2. Else, OTU non rarefied table in biom format.


27- plot_taxa
-------------
2nd step (/3) for taxonomic affiliation plot.
Make taxaonomy summary bar plots based on taxonomy assignment. 

This step takes as input files:

1. Summarized information from previous step.


28- plot_heatmap
----------------
Last step for taxonomic affiliation plot.
Make heatmap at phylum level. 

This step takes as input files:

1. Summarized information from previous step.


29- krona
---------
Plot Krona chart for taxonomic affiliation

30- plot_to_alpha
-----------------
Final report 1st part for the Amplicon-Seq pipeline. Display results (taxonomy, heatmap and alpha diversity).

31- beta_diversity
------------------
1st step (/3) for 2D PCoA plot.
Calculate beta diversity (pairwise sample dissimilarity) on OTU table. The OTU table has to be normalized. 
Only works with >= 4 samples

This step takes as input files:

1. OTU rarefied table in biom format.
2. Tree file.


32- pcoa
--------
2nd step (/3) for 2D PCoA plot.
Compute coordinates pour PCoA 

This step takes as input file:

1. Matrix produced in the previous step.


33- pcoa_plot
-------------
Last step for 2D PCoA plot.

This step takes as input file:

1. PCoA from the previous step.


34- plot_to_beta
----------------
Final report's 2nd part for the Amplicon-Seq pipeline. Display results (beta diversity PCoA plots).




Tutorial:
---------

************************************************************************
**SETTING**
************************************************************************

1a) IF you have a map file: copy the path to the "map_file" variable in the configuration file.
   eg: map_file=/path/to/map.txt  

   ELSE: Leave the variable empty. 
   eg: map_file=

1b) For the map file: All '_' character have to be replaced by '.' for the sample name. Exemple:

In readset.tsv:

Sample name : Ya_4_w3

In map file:

Sample name: Ya.4.w3

**For the following settings, you need to download the databases. Run the .sh files in resources/genomes.**

2a) For 16s study, use the Greengenes database for OTU picking. 

- $MUGQIC_INSTALL_HOME/resources/genomes/greengenes.sh

In the [DEFAULT] section:

name=greengenes

version=138

similarity_treshold=97 (or 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 99)

Also, don't forget the [qiime_otu_picking] section for the similarity treshold!

2b) For 18s study, use the Silva database.

- $MUGQIC_INSTALL_HOME/resources/genomes/silva.sh

In the [DEFAULT] section:

name=silva

version=111

similarity_treshold=97 (or 90, 94, 99)

Also, don't forget the [qiime_otu_picking] section for the similarity treshold!

2c) For ITS study, use the UNITE database for OTU picking.

- $MUGQIC_INSTALL_HOME/resources/genomes/unite.sh

In the [DEFAULT] section:

name=unite

version=1211

similarity_treshold=97 (or 99)

Also, don't forget the [qiime_otu_picking] section for the similarity treshold!

3a) For 16s and 18s study, use the GOLD database for chimera detection. 

- $MUGQIC_INSTALL_HOME/resources/genomes/chimera_gold.sh

In the [uchime] section:

name=gold
version=20110519

3b) For ITS study, use the UNITE database for chimera detection.

- Create a directory: mkdir -p $MUGQIC_INSTALL_HOME/resources/genomes/chimera_unite_db/20150311/

- You need to download it maually: https://unite.ut.ee/sh_files/uchime_reference_dataset_11_03_2015.zip

- Unzip the file and copy the database: cp uchime_sh_refs_dynamic_original_985_11.03.2015.fasta $MUGQIC_INSTALL_HOME/resources/genomes/chimera_unite_db/20150311/unite.fasta

In the [uchime] section:

name=unite

version=20150311


************************************************************************
**USAGE** 2 mods: de novo or close reference.
************************************************************************

<!-- ###################################################### -->
**de novo mod**
<!-- ###################################################### -->

1) To step 17: qiime_report

ampliconseq.py -r readsets.tsv -s 1-8,10-17 -c ampliconseq.base.ini -o analysis > analysis.sh

NB1: This step can be long, some steps can be parallelized (step 7,10,12,14 - ie: uchime, otu_picking, otu_assigning, otu_alignment)

--> Modify manually the 'ppn' to set it as same as 'threads' variable number. 

NB2: For large data, set to >= 20 threads (and ppn) for uchime.

NB3: For large data, set to >= 20 threads (and more for ppn) for otu picking.

NB4: For large data, set to >= 15 threads (and ppn) for otu_alignment.

NB5: For large data, increase the walltime (> 24h) of phylogeny step (16) 

A) Report 1

ampliconseq.py -r readset.tsv -s 1-17 -c ampliconseq.base.ini -o analysis --report > report.sh

The otu_table.sum file helps you to choose the maximum rarefaction treshold. It CAN'T be > than the max Counts/sample.

cat analysis/otus/otu_table.sum

Look for Max in Counts/sample summary:

Report to the "multiple_rarefaction_max" variable in the configuration file.

4) To step 22: qiime_report2

ampliconseq.py -r readset.tsv -s 18-22 -c ampliconseq.base.ini -o analysis > analysis2.sh

B) Report 2

ampliconseq.py -r readset.tsv -s 1-22 -c ampliconseq.base.ini -o analysis --report > report.sh

**Rarefaction normalization**

The report helps you to choose the single rarefaction treshold (for all the samples) by looking the rarefaction curves for each sample.
Report to the "single_rarefaction_depth" variable in the configuration file.

*if you have < 4 samples*

5) To step 30: plot_to_alpha

ampliconseq.py -r readset.tsv -s 23,25-30 -c ampliconseq.base.ini -o analysis > analysis3.sh

C) Report 3: analysis done

ampliconseq.py -r readset.tsv -s 1-30 -c ampliconseq.base.ini -o analysis --report > report.sh

*else*

5) To step 34: plot_to_beta

ampliconseq.py -r readset.tsv -s 23,25-34 -c ampliconseq.base.ini -o analysis > analysis3.sh

C) Report 3: analysis done

ampliconseq.py -r readset.tsv -s 1-23,25-34 -c ampliconseq.base.ini -o analysis --report > report.sh

**CSS normalization**

*if you have < 4 samples*

5) To step 30: plot_to_alpha

ampliconseq.py -r readset.tsv -s 24-30 -c ampliconseq.base.ini -o analysis > analysis3.sh

C) Report 3: analysis done

ampliconseq.py -r readset.tsv -s 1-22,24-30 -c ampliconseq.base.ini -o analysis --report > report.sh

*else*

5) To step 34: plot_to_beta

ampliconseq.py -r readset.tsv -s 24-34 -c ampliconseq.base.ini -o analysis > analysis3.sh

C) Report 3: analysis done

ampliconseq.py -r readset.tsv -s 1-22,24-34 -c ampliconseq.base.ini -o analysis --report > report.sh

<!-- ###################################################### -->
**Close ref mod (for large dataset)**
<!-- ###################################################### -->

1) To step 17: qiime_report

ampliconseq.py -r readset.tsv -s 1-9,11-17 -c ampliconseq.base.ini -o analysis > analysis.sh

NB1: This step can be long, some steps can be parallelized (step 7,9,12,14 - ie: uchime, otu_ref_picking, otu_assigning, otu_alignment)

NB2: For large data, set to >= 20 threads (and ppn) for uchime.

NB3: For large data, set to >= 20 threads (and more for ppn) for otu_ref_picking.

NB4: For large data, set to >= 15 threads (and ppn) for otu_alignment.

NB5: For large data, increase the walltime (> 24h) of phylogeny step (16) 

A) Report 1

ampliconseq.py -r readset.tsv -s 1-17 -c ampliconseq.base.ini -o analysis --report > report.sh

The otu_table.sum file helps you to choose the maximum rarefaction treshold. It CAN'T be > than the max Counts/sample.

cat analysis/otus/otu_table.sum

Look for Max in Counts/sample summary:

Report to the "multiple_rarefaction_max" variable in the configuration file.

4) To step 22: qiime_report2

ampliconseq.py -r readset.tsv -s 18-22 -c ampliconseq.base.ini -o analysis > analysis2.sh

B) Report 2

ampliconseq.py -r readset.tsv -s 1-22 -c ampliconseq.base.ini -o analysis --report > report.sh

**Rarefaction normalization**

The report helps you to choose the single rarefaction treshold (for all the samples) by looking the rarefaction curves for each sample.
Report to the "single_rarefaction_depth" variable in the configuration file.

*if you have < 4 samples*

5) To step 30: plot_to_alpha

ampliconseq.py -r readset.tsv -s 23,25-30 -c ampliconseq.base.ini -o analysis > analysis3.sh

C) Report 3: analysis done

ampliconseq.py -r readset.tsv -s 1-30 -c ampliconseq.base.ini -o analysis --report > report.sh

*else*

5) To step 34: plot_to_beta

ampliconseq.py -r readset.tsv -s 23,25-34 -c ampliconseq.base.ini -o analysis > analysis3.sh

C) Report 3: analysis done

ampliconseq.py -r readset.tsv -s 1-23,25-34 -c ampliconseq.base.ini -o analysis --report > report.sh

**CSS normalization**

*if you have < 4 samples*

5) To step 30: plot_to_alpha

ampliconseq.py -r readset.tsv -s 24-30 -c ampliconseq.base.ini -o analysis > analysis3.sh

C) Report 3: analysis done

ampliconseq.py -r readset.tsv -s 1-22,24-30 -c ampliconseq.base.ini -o analysis --report > report.sh

*else*

5) To step 34: plot_to_beta

ampliconseq.py -r readset.tsv -s 24-34 -c ampliconseq.base.ini -o analysis > analysis3.sh

C) Report 3: analysis done

ampliconseq.py -r readset.tsv -s 1-22,24-34 -c ampliconseq.base.ini -o analysis --report > report.sh


####################################################################################################################################
####################################################################################################################################
