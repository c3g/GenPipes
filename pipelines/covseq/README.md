[TOC]


CoVSeq Pipeline
================

pwet


Usage
-----
```
#!text

usage: covseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                 [--no-json] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [--sanity-check]
                 [--container {wrapper, singularity} <IMAGE PATH>]
                 [--genpipes_file GENPIPES_FILE] [-r READSETS] [-v]

Version: 4.1.2

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
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
```
![covseq workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_covseq.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_covseq.png)
```
------
1- host_reads_removal
2- kraken_analysis
3- cutadapt
4- mapping_bwa_mem_sambamba
5- sambamba_merge_sam_files
6- sambamba_filtering
7- ivar_trim_primers
8- covseq_metrics
9- freebayes_calling
10- ivar_calling
11- snpeff_annotate
12- ivar_create_consensus
13- bcftools_create_consensus
14- quast_consensus_metrics
15- rename_consensus_header_ivar
16- rename_consensus_header_freebayes
17- ncovtools_quickalign
18- prepare_table
19- prepare_report_ivar
20- prepare_report_freebayes

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
kraken

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
BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

This step takes as input files:

1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
2. Else, BAM files from the readset file

sambamba_filtering
------------------
Filter raw bams with sambamba and an awk cmd to filter by insert size

ivar_trim_primers
-----------------
Remove primer sequences to individual bam files using ivar

covseq_metrics
--------------
Multiple metrics from dnaseq

freebayes_calling
-----------------
FreeBayes is a haplotype-based variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

This method avoids one of the core problems with alignment-based variant detection that identical sequences may have multiple possible alignments.

ivar_calling
------------
ivar calling

snpeff_annotate
---------------
Consensus annotation with SnpEff

ivar_create_consensus
---------------------
Create consensus with ivar through a samtools mpileup

bcftools_create_consensus
-------------------------
bcftools consensus creation

quast_consensus_metrics
-----------------------
Generate QUAST metrics on consensus

rename_consensus_header_ivar
----------------------------
Rename reads headers

rename_consensus_header_freebayes
---------------------------------
Rename reads headers

ncovtools_quickalign
--------------------
Uses quickalign to provides summary statistics, which can be used to determine the sequencing quality and evolutionary novelty of input genomes (e.g. number of new mutations and indels). 

It uses ivar consensus as well as freebayes consensus to arrive at the alignment decisions.

prepare_table
-------------
Gathers all analysis data for quast, kraken and other metrics and module details.

prepare_report_ivar
-------------------
Prepare ivar analysis report.

prepare_report_freebayes
------------------------
Prepare ivar analysis report.


