[TOC]

ampliconseq Pipeline

Usage
-----


```
usage: genpipes ampliconseq [-h] -c CONFIG [CONFIG ...] [-s STEPS]
                            [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                            [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                            [--no-json] [--json-pt] [--clean]
                            [--container {wrapper, singularity} <IMAGE PATH>]
                            [--genpipes_file GENPIPES_FILE]
                            [-l {debug,info,warning,error,critical}]
                            [--sanity-check] [--wrap [WRAP]] -r READSETS_FILE
                            [-d DESIGN_FILE] [-v]

Version: 5.0.0

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

options:
  -h, --help            show this help message and exit
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
  --force_mem_per_cpu FORCE_MEM_PER_CPU
                        Take the mem input in the ini file and force to have a
                        minimum of mem_per_cpu by correcting the number of cpu
                        (default: None)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a validsingularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
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
0 trimmomatic16S
1 merge_trimmomatic_stats16S
2 flash_pass1
3 ampliconLengthParser
4 flash_pass2
5 merge_flash_stats
6 asva
7 multiqc
```
trimmomatic16S 
--------------
 
MiSeq raw reads adapter & primers trimming and basic QC is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', parameter 'adapter_fasta'),
it is used first. Else, Adapter1, Adapter2, Primer1 and Primer2 columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. Sequences are reversed-complemented and swapped.

This step takes as input files:
1. MiSeq paired-End FASTQ files from the readset file

merge_trimmomatic_stats16S 
--------------------------
 
The trim statistics per readset are merged at this step.

flash_pass1 
-----------
 
Merges paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/). Overlapping regions between paired-end reads are found and 
then merged into a continuous strand.

ampliconLengthParser 
--------------------
 
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

