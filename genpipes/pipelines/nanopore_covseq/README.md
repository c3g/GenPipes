[TOC]


Nanopore CoVSeq Pipeline
==============

For information on the structure and contents of the Nanopore readset file, please consult [here](https://bitbucket.org/mugqic/genpipes/src/master/#markdown-header-nanopore).
    

Usage
-----


```
#!text

usage: genpipes nanopore_covseq [-h] [--clean] -c CONFIG [CONFIG ...]
                                [--container {wrapper, singularity} <IMAGE PATH>]
                                [-f] [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                                [--genpipes_file GENPIPES_FILE]
                                [-j {pbs,batch,daemon,slurm}] [--json-pt]
                                [-l {debug,info,warning,error,critical}]
                                [--no-json] [-o OUTPUT_DIR] [--sanity-check]
                                [-s STEPS] [--wrap [WRAP]] -r READSETS_FILE
                                [-d DESIGN_FILE] [-v]
                                [-t {default,basecalling}]

Version: 5.0.0

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
  -t {default,basecalling}, --type {default,basecalling}
                        Type of CoVSeQ analysis,basecalling on/off (default
                        without basecalling)

Steps:

Protocol default
1 host_reads_removal
2 kraken_analysis
3 artic_nanopolish
4 wub_metrics
5 covseq_metrics
6 snpeff_annotate
7 quast_consensus_metrics
8 rename_consensus_header
9 prepare_report

Protocol basecalling
1 guppy_basecall
2 guppy_demultiplex
3 pycoqc
4 host_reads_removal_dependency
5 kraken_analysis
6 artic_nanopolish
7 wub_metrics
8 covseq_metrics
9 snpeff_annotate
10 quast_consensus_metrics
11 rename_consensus_header
12 prepare_report
```

host_reads_removal 
------------------
 
Runs [minimap2](https://github.com/lh3/minimap2) on a hybrid genome to remove potential host reads.

kraken_analysis 
---------------
 
Taxonomic sequence classification system using [kraken](https://github.com/DerrickWood/kraken2).

artic_nanopolish 
----------------
 
Runs artic nanopolish pipeline on all samples.

wub_metrics 
-----------
 
Generate WUB metrics on bam file.

covseq_metrics 
--------------
 
Collect metrics.

snpeff_annotate 
---------------
 
Consensus annotation with [SnpEff](https://pcingola.github.io/SnpEff/).

quast_consensus_metrics 
-----------------------
 
Generate [QUAST](http://quast.sourceforge.net/) metrics on consensus.

rename_consensus_header 
-----------------------
 
Rename reads headers.

prepare_report 
--------------
 
Prepare analysis report.

guppy_basecall 
--------------
 
Use guppy to perform basecalling on raw FAST5 files.

guppy_demultiplex 
-----------------
 
Use guppy to perform demultiplexing on raw FASTQ read files.

pycoqc 
------
 
Use [pycoQC](https://hpc.nih.gov/apps/pycoQC.html) to produce an interactive quality report based on the summary file and
alignment outputs.

host_reads_removal_dependency 
-----------------------------
 
Runs [minimap2](https://github.com/lh3/minimap2) on a hybrid genome to remove potential host reads.

