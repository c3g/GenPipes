[TOC]



Usage
-----
```
#!text

usage: rnaseq_light.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                       [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                       [--no-json] [--report] [--clean]
                       [-l {debug,info,warning,error,critical}]
                       [--sanity-check]
                       [--container {wrapper, singularity} <IMAGE PATH>]
                       [--genpipes_file GENPIPES_FILE] [-d DESIGN]
                       [-r READSETS] [-v]

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
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
```
![rnaseq_light workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq_light.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq_light.png)
```
------
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- kallisto
5- kallisto_count_matrix
6- gq_seq_utils_exploratory_analysis_rnaseq_light
7- sleuth_differential_expression

```
picard_sam_to_fastq
-------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

trimmomatic
-----------
Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
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
Run Kallisto on fastq files for a fast esimate of abundance.

kallisto_count_matrix
---------------------
Use the output from Kallisto to create a transcript count matrix. 

gq_seq_utils_exploratory_analysis_rnaseq_light
----------------------------------------------
Exploratory analysis using the gqSeqUtils R package adapted for RnaSeqLight

sleuth_differential_expression
------------------------------
Performs differential gene expression analysis using [Sleuth](http://pachterlab.github.io/sleuth/). 
Analysis are performed both at a transcript and gene level, using two different tests: LRT and WT. 


