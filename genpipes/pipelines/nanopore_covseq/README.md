[TOC]


    Nanopore CoVSeq Pipeline
    ==============

    For information on the structure and contents of the Nanopore readset file, please consult [here](
    https://bitbucket.org/mugqic/genpipes/src/master/#markdown-header-nanopore).
    
================

Usage
-----


```
#!text

usage: genpipes nanopore_covseq [-h] -c CONFIG [CONFIG ...] [-s STEPS]
                                [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}]
                                [-f] [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                                [--no-json] [--json-pt] [--clean]
                                [--container {wrapper, singularity} <IMAGE PATH>]
                                [--genpipes_file GENPIPES_FILE]
                                [-l {debug,info,warning,error,critical}]
                                [--sanity-check] [--wrap [WRAP]] -r
                                READSETS_FILE [-d DESIGN_FILE] [-v]
                                [-t {default,basecalling}]

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
  -t {default,basecalling}, --type {default,basecalling}
                        Type of CoVSeQ analysis,basecalling on/off (default
                        without basecalling)

Steps:

Protocol default
0 host_reads_removal
1 kraken_analysis
2 artic_nanopolish
3 wub_metrics
4 covseq_metrics
5 snpeff_annotate
6 quast_consensus_metrics
7 rename_consensus_header
8 prepare_report

Protocol basecalling
0 guppy_basecall
1 guppy_demultiplex
2 pycoqc
3 host_reads_removal_dependency
4 kraken_analysis
5 artic_nanopolish
6 wub_metrics
7 covseq_metrics
8 snpeff_annotate
9 quast_consensus_metrics
10 rename_consensus_header
11 prepare_report
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

