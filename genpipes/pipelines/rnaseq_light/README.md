[TOC]

RNA-Seq Light Pipeline
================

The RNA-Seq Light Pipeline is a lightweight analysis of gene expression in RNA sequencing data. 
The pipeline is based on [Kallisto](https://pachterlab.github.io/kallisto/about.html) and differential expression analysis is performed by [Sleuth](http://pachterlab.github.io/sleuth/).
It is especially useful for quick Quality Control (QC) in gene sequencing studies.
    
Usage
-----

```
#!text
usage: genpipes rnaseq_light [-h] [--clean] -c CONFIG [CONFIG ...]
                             [--container {wrapper, singularity} <IMAGE PATH>]
                             [-f] [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                             [--genpipes_file GENPIPES_FILE]
                             [-j {pbs,batch,daemon,slurm}] [--json-pt]
                             [-l {debug,info,warning,error,critical}]
                             [--no-json] [-o OUTPUT_DIR] [--sanity-check]
                             [-s STEPS] [--wrap [WRAP]] -r READSETS_FILE
                             [-d DESIGN_FILE] [-v]

Version: 5.0.3

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

options:
  -h, --help            show this help message and exit
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
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
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -j {pbs,batch,daemon,slurm}, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  --wrap [WRAP]         Path to the genpipe cvmfs wrapper script. Default is g
                        enpipes/ressources/container/bin/container_wrapper.sh.
                        This is a convenience options for using genpipes in a
                        container
  -r READSETS_FILE, --readsets READSETS_FILE
                        readset file
  -d DESIGN_FILE, --design DESIGN_FILE
                        design file
  -v, --version         show the version information and exit

Steps:

Protocol default
1 picard_sam_to_fastq
2 trimmomatic
3 merge_trimmomatic_stats
4 kallisto
5 kallisto_count_matrix
6 gq_seq_utils_exploratory_analysis_rnaseq_light
7 sleuth_differential_expression
8 multiqc
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

kallisto 
--------
 
Run [Kallisto](https://pachterlab.github.io/kallisto/about.html) on fastq files for a fast esimate of abundance.

kallisto_count_matrix 
---------------------
 
Use the output from Kallisto to create a transcript count matrix.
Create a summary table to be included in the multiqc report.

gq_seq_utils_exploratory_analysis_rnaseq_light 
----------------------------------------------
 
Exploratory analysis using the gqSeqUtils R package adapted for RnaSeqLight.

sleuth_differential_expression 
------------------------------
 
Performs differential gene expression analysis using [Sleuth](http://pachterlab.github.io/sleuth/).
Analysis are performed both at a transcript and gene level, using two different tests: LRT and WT.

multiqc 
-------
 
A quality control report for all samples is generated.
For more detailed information about MultiQC visit: [MultiQC](http://multiqc.info/)

