[TOC]


CoVSeQ Nanopore Pipeline
==============
The SOP for Nanopore data is based on the ARTIC SARS-CoV2 pipeline using nanopolish. Their full documentation is found [here](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

The protocol was closely followed with the majority of changes, involving technical adaptations to be able to run in a High Performance Computing environment where the usage Conda is not advisable. In summary, if basecalling protocol is selected, the pipeline will do basecalling using Guppy (GPU) and demultiplexing. Then, for all samples, the pipeline will do de-hosting, run the `artic-nanopolish` wrapper which performs alignment to the SARS-CoV2 reference (using `minimap2`), variant calling (using `nanopolish`), and consensus generation (using `artic_mask` + `bcftools consensus`). Finally, custom scripts and `ncov_tools` are run to report on quality metrics. 

*Important note*: the pipeline is set up to use ARTIC v3 amplicon scheme as a default. If ARTIC v4 is required, use the appropriate `.ini` file. For all other amplicon schemes, add the appropriate primer and amplicon bed files and use a custom `.ini` for processing. 


Usage
-----
```
#!text

usage: nanopore_covseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                   [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                   [--no-json] [--report] [--clean]
                   [-l {debug,info,warning,error,critical}] [--sanity-check]
                   [-t {default,basecalling}]
                   [--container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}]
                   [-r READSETS] [-v]

Version: 4.0.0

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
                        job scheduler type (default: slurm)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  -t {stringtie,cufflinks}, --type {stringtie,cufflinks} CoVSeQ analysis type
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  --container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}
                        run pipeline inside a container providing a container
                        image path or accessible docker/singularity hub path
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

---
default:
1- host_reads_removal,
2- kraken_analysis,
3- artic_nanopolish,
4- wub_metrics,
5- covseq_metrics,
6- snpeff_annotate,
7- quast_consensus_metrics,
8- rename_consensus_header,
9- prepare_report

---
basecalling:
1- guppy_basecall,
2- guppy_demultiplex,
3- pycoqc,
4- host_reads_removal,
5- kraken_analysis,
6- artic_nanopolish,
7- wub_metrics,
8- covseq_metrics,
9- snpeff_annotate,
10- quast_consensus_metrics,
11- rename_consensus_header,
12- prepare_report

```

guppy_basecall
-------
Use the Oxford Nanopore basecaller Guppy to basecall raw FAST5 files and produce FASTQ files. Basecalling model `dna_r9.4.1_450bps_hac.cfg` is used by default. 

guppy_demultiplex
-------
Use the Ofxord Nanopore basecaller Guppy to demultiplex FASTQ files based on their barcode. Barcode arrangement `barcode_arrs_nb96.cfg` is used by default. **Important** the parameter `--require_barcodes_both_ends` is set by default. 

pycoqc
-------
If basecalling and demultiplexing were performed, a `pycoQC` interactive report is produced to aid with the run QC. 

host_reads_removal
-------
Using a mapping approach with a hybrid GRCh38 + SARS-CoV2 genome, reads that map to the Human Genome are removed from the analysis. A "de-hosted" FASTQ is produced. 

kraken_analysis
-------
Additionally, `kraken2` is used to produce a report on the raw data, which can be used to detect additional host contamination. 

artic_nanopolish
-------
The `artic nanopolish` pipeline is used to produce consensus sequences and VCFs. Since `nanopolish` is used, this step requires both FAST5 and FASTQ files.

wub_metrics
-------
Alignment metrics are calculated using the tool `wub`. 

covseq_metrics
-------
Using all previous metrics calculated so far, a table is produced with a summary of all metrics for each individual sample. 

snpeff_annotate
-------
The VCF produced by `artic_nanopolish` is annotated using `SnpEff`. 

quast_consensus_metrics
-------
Consensus metrics are calculated using the tool `QUAST`

rename_consensus_header
-------
A final consensus sequence is produced, with the appropriate header and naming convention based on genome completeness. 

prepare_report
-------
Using `ncov_tools` and additional R scripts, final reports are produced for all samples in the run, including basic QC plots as well as a preliminary lineage assignment (as a part of `ncov_tools`). 
