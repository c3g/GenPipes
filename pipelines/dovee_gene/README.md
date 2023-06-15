[TOC]


DOvEE Gene Pipeline


Usage
-----
```
#!text

usage: dovee_gene.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                     [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                     [--no-json] [--report] [--clean]
                     [-l {debug,info,warning,error,critical}] [--sanity-check]
                     [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                     [--container {wrapper, singularity} <IMAGE PATH>]
                     [--genpipes_file GENPIPES_FILE] [-d DESIGN] [-p PAIRS]
                     [-t {vardict,copy-number}] [-r READSETS] [-v]

Version: 4.3.2

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
  --force_mem_per_cpu FORCE_MEM_PER_CPU
                        Take the mem input in the ini file and force to have a
                        minimum of mem_per_cpu by correcting the number of cpu
                        (default: None)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -d DESIGN, --design DESIGN
                        File indicating whether sample is saliva or brush
  -p PAIRS, --pairs PAIRS
                        File with sample pairing information
  -t {vardict,copy-number}, --type {vardict,copy-number}
                        Type of pipeline (default vardict)
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
vardict:
1- trimmer
2- fastp
3- bwa_mem_samtools_sort
4- samtools_merge
5- locatit_dedup_bam
6- samtools_subsample
7- mosdepth
8- picard_metrics
9- vardict_single
10- bcftools_stats
11- conpair_concordance
12- multiqc
----
copy-number:
1- trimmer
2- fastp
3- bwa_mem_samtools_sort
4- samtools_merge
5- locatit_hybrid_dedup
6- mosdepth
7- hmm_readCounter
8- run_ichorCNA
9- multiqc

```

trimmer
-------
Read trimming with AGeNT Trimmer from Agilent Genomics NextGen tool kit. 

fastp
-----
Generate basic QC metrics for trimmed reads with fastp.

bwa_mem_samtools_sort
---------------------
The trimmed reads are aligned to the reference genome with bwa mem, followed by sorting with samtools.
The bwa mem output is piped directory into samtools to avoid saving the intermediate file.

samtools_merge
--------------
Merges sorted bam files for each sample, if there are multiple readsets for the sample. If only a single readset file is found, a symlink to the bam is created. 

locatit_dedup_bam
-----------------
Deduplicate bam files with AGeNT locatit in hybrid (saliva only) or duplex (both saliva and brush) mode.
Used for SureSelect samples in vardict protocol.
Saliva and brush samples are identified via the design file. 

samtools_subsample
------------------
Vardict protocol only.
Subsample the deduplicated bams to 1.5M reads if the number of overlapping reads exceeds this number (brush samples only).
Number of overlapping reads is assessed with samtools view -c, from which a fraction is calculated for subsampling.
If fraction >= 1, the bam is not subsampled.

mosdepth
--------
Calculate depth stats for captured regions with mosdepth.

picard_metrics
--------------
Collect on and off target metrics for SureSelect samples with Picard.
Collect gc bias metrics for SureSelect samples with Picard.

vardict_single
--------------
Variant calling with vardict in single mode.

bcftools_stats
--------------
Collect stats number and types of variants in vcf files produced by vardict_single.

conpair_concordance
-------------------
Conpair is a fast and robust method dedicated for human tumor-normal studies to perform concordance verification
(= samples coming from the same individual).
Run once per patient to ensure that saliva and brush samples are correctly assigned.
Requires pair file.

multiqc
-------
Aggregate results from bioinformatics analyses across many samples into a single report
MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).

locatit_hybrid_dedup
--------------------
Deduplicate bam files with locatit in hybrid mode only. Used for low pass copy-number protocol.

hmm_readCounter
---------------
Counting number of reads in non-overlapping windows of fixed width directly from BAM files with HMM Copy Utils readCounter:
https://github.com/shahcompbio/hmmcopy_utils
Run on both saliva and brush, as both wigs needed as input for IchorCNA.

run_ichorCNA
------------
https://github.com/broadinstitute/ichorCNA
'Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.'
Run once per patient using paired saliva and brush sample wigs, requires knowing which saliva and brush came from same patient.


