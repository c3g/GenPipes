[TOC]


Nanopore CoVSeq Pipeline
==============


For information on the structure and contents of the Nanopore readset file, please consult [here](
https://bitbucket.org/mugqic/genpipes/src/master/#markdown-header-nanopore).


Usage
-----
```
#!text

usage: nanopore_covseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                          [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                          [--no-json] [--report] [--clean]
                          [-l {debug,info,warning,error,critical}]
                          [--sanity-check]
                          [--container {wrapper, singularity} <IMAGE PATH>]
                          [--genpipes_file GENPIPES_FILE] [-r READSETS]
                          [-t {default,basecalling}] [-v]

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
  -t {default,basecalling}, --type {default,basecalling}
                        Type of CoVSeQ analysis,basecalling on/off (default
                        without basecalling)
  -v, --version         show the version information and exit

Steps:
------

----
```
![nanopore_covseq default workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_nanopore_covseq_default.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_nanopore_covseq_default.png)
```
default:
1- host_reads_removal
2- kraken_analysis
3- artic_nanopolish
4- wub_metrics
5- covseq_metrics
6- snpeff_annotate
7- quast_consensus_metrics
8- rename_consensus_header
9- prepare_report
----
```
![nanopore_covseq basecalling workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_nanopore_covseq_basecalling.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_nanopore_covseq_basecalling.png)
```
basecalling:
1- guppy_basecall
2- guppy_demultiplex
3- pycoqc
4- host_reads_removal_dependency
5- kraken_analysis
6- artic_nanopolish
7- wub_metrics
8- covseq_metrics
9- snpeff_annotate
10- quast_consensus_metrics
11- rename_consensus_header
12- prepare_report

```
host_reads_removal
------------------
Runs minimap2 on a hybrid genome to remove potential host reads

kraken_analysis
---------------
kraken

artic_nanopolish
----------------
Runs artic nanopolish pipeline on all samples.

wub_metrics
-----------
Generate WUB metrics on bam file

covseq_metrics
--------------

snpeff_annotate
---------------
Consensus annotation with SnpEff

quast_consensus_metrics
-----------------------
Generate QUAST metrics on consensus

rename_consensus_header
-----------------------
Rename reads headers

prepare_report
--------------
guppy_basecall
--------------
Use guppy to perform basecalling on raw FAST5 files


guppy_demultiplex
-----------------
Use guppy to perform demultiplexing on raw FASTQ read files


pycoqc
------
Use pycoQC to produce an interactive quality report based on the summary file and
alignment outputs.

host_reads_removal_dependency
-----------------------------
Runs minimap2 on a hybrid genome to remove potential host reads


