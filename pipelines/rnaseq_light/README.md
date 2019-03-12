[TOC]


RNA-Seq Light Pipeline
======================

This alternative RNA-Seq pipeline is designed for fast, easy estimation of transcript abundances using the [Kallisto](http://pachterlab.github.io/kallisto) pseudoaligner. It requires a fully annotated refrence genome and will only estimate abundances of transcripts that are fully annotated. The [Sleuth](http://pachterlab.github.io/sleuth/) R-package is used perform differential transcript and gene analysis. The main use of this pipeline is for simple experiments on model organisms, where *fast turnaround time* and *efficient use of computing resources* are important. It is _not_ recommended for complex analyses, since it does not produce alignment files, nor does it support novel transcript detection.   

Currently, this pipeline only fully supports *paired-end libraries* with only *one readset per sample* and at least *3 biological replicates* per exeprimental group. Additional functionality will be added in future versions. 

**Note**: As of version 3.1.4, only Homo_sapiens.GRCh37 and Mus_musculus.GRCm38 are supported by default. Support for additional genomes will be added gradually. Users can also generate the required tx2gene files for their genome of preference to use this pipeline.  

Usage
-----
```
#!text

usage: rnaseq_light.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                       [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                       [--json] [--report] [--clean]
                       [-l {debug,info,warning,error,critical}] [-d DESIGN]
                       [-r READSETS] [-v]

Version: 3.1.3

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
  -j {pbs,batch,daemon,slurm}, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: pbs)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --json                create a JSON file per analysed sample to track the
                        analysis status (default: false)
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
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
```
```
------
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- kallisto
5- kallisto_count_matrix
6- gq_seq_utils_exploratory_analysis_rnaseq_light
7- sleuth_differential_expresssion

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
Run [Kallisto](http://pachterlab.github.io/kallisto) on fastq files for a fast esimate of abundance.

kallisto_count_matrix
---------------------
Generate a kallisto pseudo-count matrix with appropriate format for following steps.

gq_seq_utils_exploratory_analysis_rnaseq_light
----------------------------------------------
Exploratory analysis using the gqSeqUtils R package adapted for RnaSeqLight

sleuth_differential_expression
------------------------------
Use [Sleuth](http://pachterlab.github.io/sleuth/) to perform differential gene and transcript expression. Requires a reference tx2gene file. 

