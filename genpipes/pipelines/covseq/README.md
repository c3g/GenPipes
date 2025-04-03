<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [CoVSeq Pipeline](#covseq-pipeline)
  - [The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage](#the-pipeline-is-designed-to-be-run-on-a-cluster-and-is-configured-using-a-configuration-file-the-pipeline-can-be-run-in-a-single-step-or-in-multiple-steps-the-pipeline-can-also-be-run-in-parallel-to-process-multiple-samples-simultaneously%0Ausage)
  - [host_reads_removal](#host_reads_removal)
  - [kraken_analysis](#kraken_analysis)
  - [cutadapt](#cutadapt)
  - [mapping_bwa_mem_sambamba](#mapping_bwa_mem_sambamba)
  - [sambamba_merge_sam_files](#sambamba_merge_sam_files)
  - [sambamba_filtering](#sambamba_filtering)
  - [ivar_trim_primers](#ivar_trim_primers)
  - [covseq_metrics](#covseq_metrics)
  - [freebayes_calling](#freebayes_calling)
  - [ivar_calling](#ivar_calling)
  - [snpeff_annotate](#snpeff_annotate)
  - [ivar_create_consensus](#ivar_create_consensus)
  - [bcftools_create_consensus](#bcftools_create_consensus)
  - [quast_consensus_metrics](#quast_consensus_metrics)
  - [rename_consensus_header_ivar](#rename_consensus_header_ivar)
  - [rename_consensus_header_freebayes](#rename_consensus_header_freebayes)
  - [ncovtools_quickalign](#ncovtools_quickalign)
  - [prepare_table](#prepare_table)
  - [prepare_report_ivar](#prepare_report_ivar)
  - [prepare_report_freebayes](#prepare_report_freebayes)
  - [run_multiqc](#run_multiqc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

CoVSeq Pipeline
================
A pipeline to process and analyze SARS-CoV-2 sequencing data from Illumina platforms. The pipeline uses Cutadapt for adapter trimming, Kraken for taxonomic classification, BWA for read alignment, Sambamba for sorting and indexing, and Freebayes for variant calling. The pipeline also includes a number of metrics to assess the quality.

The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage
-----

```
#!text
usage: genpipes covseq [-h] [--clean] -c CONFIG [CONFIG ...]
                       [--container {wrapper, singularity} <IMAGE PATH>] [-f]
                       [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                       [--genpipes_file GENPIPES_FILE]
                       [-j {pbs,batch,daemon,slurm}] [--json-pt]
                       [-l {debug,info,warning,error,critical}]
                       [-o OUTPUT_DIR] [--sanity-check] [-s STEPS]
                       [--wrap [WRAP]] -r READSETS_FILE [-d DESIGN_FILE] [-v]

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

Steps:

Protocol default
1 host_reads_removal
2 kraken_analysis
3 cutadapt
4 mapping_bwa_mem_sambamba
5 sambamba_merge_sam_files
6 sambamba_filtering
7 ivar_trim_primers
8 covseq_metrics
9 freebayes_calling
10 ivar_calling
11 snpeff_annotate
12 ivar_create_consensus
13 bcftools_create_consensus
14 quast_consensus_metrics
15 rename_consensus_header_ivar
16 rename_consensus_header_freebayes
17 ncovtools_quickalign
18 prepare_table
19 prepare_report_ivar
20 prepare_report_freebayes
21 run_multiqc
```

host_reads_removal 
------------------
 
The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

This step takes as input files:
1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

kraken_analysis 
---------------
 
Taxonomic sequence classification system using [kraken](https://github.com/DerrickWood/kraken2).

cutadapt 
--------
 
Raw reads quality trimming and removing adapters is performed using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html).
'Adapter1' and 'Adapter2' columns from the readset file ar given to Cutadapt. For PAIRED_END readsets, both adapters are used.
For SINGLE_END readsets, only Adapter1 is used and left unchanged.
To trim the front of the read use adapter_5p_fwd and adapter_5p_rev (for PE only) in cutadapt section of ini file.

This step takes as input files:
1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

mapping_bwa_mem_sambamba 
------------------------
 
The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

This step takes as input files:
1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

sambamba_merge_sam_files 
------------------------
 
BAM readset files are merged into one file per sample. Merge is done using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

This step takes as input files:

1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
2. Else, BAM files from the readset file

sambamba_filtering 
------------------
 
Filter raw bams with [Sambamba](http://lomereiter.github.io/sambamba/index.html) and an awk cmd to filter by insert size

ivar_trim_primers 
-----------------
 
Remove primer sequences to individual bam files using [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html).

covseq_metrics 
--------------
 
Gathering multiple metrics.

freebayes_calling 
-----------------
 
[FreeBayes](https://github.com/freebayes/freebayes) is a haplotype-based variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

This method avoids one of the core problems with alignment-based variant detection that identical sequences may have multiple possible alignments.

ivar_calling 
------------
 
[ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) calling variants.

snpeff_annotate 
---------------
 
Consensus annotation with [SnpEff](https://pcingola.github.io/SnpEff/).

ivar_create_consensus 
---------------------
 
Create consensus with [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) through a [samtools](http://samtools.sourceforge.net/) mpileup

bcftools_create_consensus 
-------------------------
 
[bcftools](https://samtools.github.io/bcftools/bcftools.html) consensus creation

quast_consensus_metrics 
-----------------------
 
Generate [QUAST](http://quast.sourceforge.net/) metrics on consensus

rename_consensus_header_ivar 
----------------------------
 
Rename reads headers after [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) variants calling.

rename_consensus_header_freebayes 
---------------------------------
 
Rename reads headers after [FreeBayes](https://github.com/freebayes/freebayes) variants calling.

ncovtools_quickalign 
--------------------
 
Uses [ncov-tools](https://github.com/jts/ncov-tools) quickalign to provides summary statistics, which can be used to determine the sequencing quality and evolutionary novelty of input genomes (e.g. number of new mutations and indels). 

It uses ivar consensus as well as freebayes consensus to arrive at the alignment decisions.

prepare_table 
-------------
 
Gathers all analysis data for [QUAST](http://quast.sourceforge.net/), [kraken](https://github.com/DerrickWood/kraken2) and other metrics and module details.

prepare_report_ivar 
-------------------
 
Prepare [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) analysis report.

prepare_report_freebayes 
------------------------
 
Prepare [FreeBayes](https://github.com/freebayes/freebayes) analysis report.

run_multiqc 
-----------
 
Run [multiqc](https://multiqc.info/) on all samples.

