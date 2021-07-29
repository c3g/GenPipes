[TOC]

Salmon Proof-of-Concept RNA-seq Pipeline
================
This pipeline is a minimal example implementation for a bioinformatics RNA-seq 
workflow using GenPipes. It is based on a basic Salmon quantification analysis
as shown in [this repository](https://github.com/GoekeLab/bioinformatics-workflows).

The implementation consists of three steps, the first is a basic quality control using
fastQC. The second and third steps use Salmon to index a reference transcriptome and 
quantify transcripts in a sample. 

The implementation assumes the presence of linux modules for fastQC, Salmon and Java i
(OpenJDK). As with all GenPipes pipelines, it supports running as a batch or using 
the PBS and SLURM schedulers. 

Standard GenPipes features like dependency management, readset/sample merging, and 
job reporting are available. Since this is a minimal implementation, additional features
such such as protocols, designs, etc. are not available for this pipeline. 

For more information, please consult the [GenPipes Documentation](https://genpipes.readthedocs.io/en/master/index.html)

Usage
-----
```
#!text

usage: goeke_rnaseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                       [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                       [--no-json] [--report] [--clean]
                       [-l {debug,info,warning,error,critical}]
                       [--sanity-check]
                       [--container {wrapper, singularity} <IMAGE PATH>]
                       [--genpipes_file GENPIPES_FILE] [-d DESIGN]
                       [-t {cufflinks,stringtie}] [-r READSETS] [-v]

Version: 3.5.0

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
  -t {cufflinks,stringtie}, --type {cufflinks,stringtie}
                        Type of RNA-seq assembly method (default stringtie)
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
```
```
------
1- fastqc
2- salmon_index
3- salmon_quant

```
fastqc
-------------------
Step 1: Quality Control (with FastQC)

salmon_index
-----------
iStep 2: Index Creation (with Salmon)

salmon_quant
-----------------------
Step 3: Quantification (with Salmon)



