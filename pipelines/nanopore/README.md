[TOC]


Nanopore Pipeline
==============

The Nanopore is used to analyse long reads produced by the Oxford Nanopore Technologies (ONT) sequencers.
Currently, the pipeline uses minimap2 to align reads to the reference genome. Additionally, it produces
a QC report that includes an interactive dashboard with data from the basecalling summary file as well
as the alignment. A step aligning random reads to the NCBI nt database and reporting the species of the
highest hits is also done as QC.

Once the QC and alignments have been produced, Picard is used to merge readsets coming from the same
sample. Finally, SVIM is used to detect Structural Variants (SV) including deletions, insertions and
translocations. For a full summary of the types of SVs detected, please consult the following [site](
https://github.com/eldariont/svim#background-on-structural-variants-and-long-reads).

The SV calls produced by SVIM are saved as VCFs for each sample, which can then be used in downstream
analyses. No filtering is performed on the SV calls.

This pipeline currently does not perform base calling and requires both FASTQ and a sequencing_summary
file produced by a ONT supported basecaller (we recommend Guppy). Additionally, the testing and
development of the pipeline were focused on genomics applications, and functionality has not been tested
for transcriptomics or epigenomics datasets.

For more information on using ONT data for structural variant detection, as well as an alternative
approach, please consult [this GitHub repository](https://github.com/nanoporetech/pipeline-structural-variation).

For information on the structure and contents of the Nanopore readset file, please consult [here](
https://bitbucket.org/mugqic/genpipes/src/master/#markdown-header-nanopore).


Usage
-----
```
#!text

usage: nanopore.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                   [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                   [--no-json] [--report] [--clean]
                   [-l {debug,info,warning,error,critical}] [--sanity-check]
                   [--container {wrapper, singularity} <IMAGE PATH>]
                   [--genpipes_file GENPIPES_FILE] [-r READSETS] [-v]

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
  -v, --version         show the version information and exit

Steps:
```
![nanopore workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_nanopore.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_nanopore.png)
```
------
1- blastqc
2- minimap2_align
3- pycoqc
4- picard_merge_sam_files
5- svim

```
blastqc
-------
Uses BLAST to perform a basic QC test by aligning 1000bp of randomly selected
reads to the NCBI nt database in order to detect potential contamination.

minimap2_align
--------------
Uses minimap2 to align the Fastq reads that passed the minimum QC threshold to
the provided reference genome. By default, it aligns to GRCh38.

pycoqc
------
Use pycoQC to produce an interactive quality report based on the summary file and
alignment outputs.

picard_merge_sam_files
----------------------
BAM readset files are merged into one file per sample.
Merge is done using [Picard](http://broadinstitute.github.io/picard/).

This step takes as input files:
Aligned and sorted BAM output files from previous minimap2_align step

svim
----
Use SVIM to perform SV calling on each sample.


