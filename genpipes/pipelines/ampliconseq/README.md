<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Amplicon-Seq Pipeline](#amplicon-seq-pipeline)
  - [The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage](#the-pipeline-is-designed-to-be-run-on-a-cluster-and-is-configured-using-a-configuration-file-the-pipeline-can-be-run-in-a-single-step-or-in-multiple-steps-the-pipeline-can-also-be-run-in-parallel-to-process-multiple-samples-simultaneously%0Ausage)
  - [trimmomatic16S](#trimmomatic16s)
  - [merge_trimmomatic_stats16S](#merge_trimmomatic_stats16s)
  - [flash_pass1](#flash_pass1)
  - [amplicon_length_parser](#amplicon_length_parser)
  - [flash_pass2](#flash_pass2)
  - [merge_flash_stats](#merge_flash_stats)
  - [asva](#asva)
  - [multiqc](#multiqc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

Amplicon-Seq Pipeline
=================

A pipeline to process amplicon sequencing data. The pipeline is designed to handle both paired-end and single-end reads and can be used to process data from any Illumina sequencer. The pipeline uses Trimmomatic to trim adapters and primers, FLASh to merge paired-end reads, and DADA2 to infer sequence variants of microbial communities.

The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage
-----

```
#!text
usage: genpipes ampliconseq [-h] [--clean] -c CONFIG [CONFIG ...]
                            [--container {wrapper, singularity} <IMAGE PATH>]
                            [-f] [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                            [--genpipes_file GENPIPES_FILE]
                            [-j {pbs,batch,daemon,slurm}] [--json-pt]
                            [-l {debug,info,warning,error,critical}]
                            [-o OUTPUT_DIR] [--sanity-check] [-s STEPS]
                            [--wrap [WRAP]] -r READSETS_FILE [-d DESIGN_FILE]
                            [-v]

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
1 trimmomatic16S
2 merge_trimmomatic_stats16S
3 flash_pass1
4 amplicon_length_parser
5 flash_pass2
6 merge_flash_stats
7 asva
8 multiqc
```

trimmomatic16S 
--------------
 
MiSeq raw reads adapter & primers trimming and basic QC is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', parameter 'adapter_fasta'),
it is used first. Else, Adapter1, Adapter2, Primer1 and Primer2 columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. Sequences are reversed-complemented and swapped.

This step takes as input files MiSeq paired-End FASTQ files from the readset file.

merge_trimmomatic_stats16S 
--------------------------
 
The trim statistics per readset are merged at this step.

flash_pass1 
-----------
 
Merges paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/). Overlapping regions between paired-end reads are found and 
then merged into a continuous strand.

amplicon_length_parser 
----------------------
 
Looks at FLASH output statistics to set input amplicon lengths for dada2. Minimum lengths are set by ensuring that they represent 
at least 1% of the total number of amplicons.

flash_pass2 
-----------
 
Merges paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/). The second pass uses statistics obtained from the first pass
to adjust merging.

merge_flash_stats 
-----------------
 
Merges statistics from both flash passes.

asva 
----
 
Checks for design files required for PCA plots, sets up directories, links readset fastq files, and initiates 
[DADA2](https://benjjneb.github.io/dada2/). 

DADA2 is used to infer sequence variants of microbial communities.

multiqc 
-------
 
A quality control report for all samples is generated.
For more detailed information about MultiQC visit: [MultiQC](http://multiqc.info/)

