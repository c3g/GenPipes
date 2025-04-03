<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [LongRead DNA-Seq Pipeline](#longread-dna-seq-pipeline)
  - [Usage](#usage)
  - [blastqc](#blastqc)
  - [minimap2_align](#minimap2_align)
  - [pycoqc](#pycoqc)
  - [picard_merge_sam_files](#picard_merge_sam_files)
  - [svim](#svim)
  - [metrics_nanoplot](#metrics_nanoplot)
  - [pbmm2_align](#pbmm2_align)
  - [metrics_mosdepth](#metrics_mosdepth)
  - [set_deepvariant_regions](#set_deepvariant_regions)
  - [deepvariant](#deepvariant)
  - [merge_filter_deepvariant](#merge_filter_deepvariant)
  - [hificnv](#hificnv)
  - [trgt_genotyping](#trgt_genotyping)
  - [sawfish](#sawfish)
  - [annotSV](#annotsv)
  - [hiphase](#hiphase)
  - [report_cpsr](#report_cpsr)
  - [multiqc](#multiqc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

LongRead DNA-Seq Pipeline
==============

The LongRead Pipeline is used to analyse long reads produced by the Oxford Nanopore Technologies (ONT) 
and PacBio Revio sequencers. The protocols used are nanopore and revio, respectively.

Currently, the nanopore protocol of the pipeline uses minimap2 to align reads to the reference genome.
Additionally, it produces a QC report that includes an interactive dashboard with data from the basecalling
summary file as well as the alignment. A step aligning random reads to the NCBI nt database and reporting 
the species of the highest hits is also done as QC.

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

The Revio protocol uses pbmm2 to align reads to the reference genome, followed by variant calling with DeepVariant
and structural variant calling with HiFiCNV, TRGT, and Sawfish. Variants are annotated with AnnotSV and phased
with HiPhase. A CPSR report can be produced from the phased variants. Metrics on the raw and mapped reads are
collected with NanoPlot and mosdepth, respectively. 

Both protocols require as input a readset file, which provides sample metadata and paths to input data (FASTQ, FAST5 or BAM).
For information on the structure and contents of the LongRead readset file, please consult [here](https://genpipes.readthedocs.io/en/latest/get-started/concepts/readset_file.html).
    
Usage
-----

```
#!text
usage: genpipes longread_dnaseq [-h] [--clean] -c CONFIG [CONFIG ...]
                                [--container {wrapper, singularity} <IMAGE PATH>]
                                [-f] [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                                [--genpipes_file GENPIPES_FILE]
                                [-j {pbs,batch,daemon,slurm}] [--json-pt]
                                [-l {debug,info,warning,error,critical}]
                                [-o OUTPUT_DIR] [--sanity-check] [-s STEPS]
                                [--wrap [WRAP]] -r READSETS_FILE
                                [-d DESIGN_FILE] [-v] [-t {nanopore,revio}]

For more documentation, visit our website: https://genpipes.readthedocs.io

options:
  -h, --help            show this help message and exit
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -c, --config CONFIG [CONFIG ...]
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
  --genpipes_file, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -j, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
  -l, --log {debug,info,warning,error,critical}
                        log level (default: info)
  -o, --output-dir OUTPUT_DIR
                        output directory (default: current)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  -s, --steps STEPS     step range e.g. '1-5', '3,6,7', '2,4-8'
  --wrap [WRAP]         Path to the genpipes cvmfs wrapper script. Default is 
                        genpipes/ressources/container/bin/container_wrapper.sh
                        . This is a convenience option for using genpipes in a
                        container
  -r, --readsets READSETS_FILE
                        readset file
  -d, --design DESIGN_FILE
                        design file
  -v, --version         show the version information and exit
  -t, --type {nanopore,revio}
                        Type of pipeline (default nanopore)

Steps:

Protocol nanopore
1 blastqc
2 minimap2_align
3 pycoqc
4 picard_merge_sam_files
5 svim

Protocol revio
1 metrics_nanoplot
2 pbmm2_align
3 picard_merge_sam_files
4 metrics_mosdepth
5 set_deepvariant_regions
6 deepvariant
7 merge_filter_deepvariant
8 hificnv
9 trgt_genotyping
10 sawfish
11 annotSV
12 hiphase
13 report_cpsr
14 multiqc
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
Aligned and sorted BAM output files from previous minimap2_align or pbmm2_align step

svim 
----
 
Use SVIM to perform SV calling on each sample.

metrics_nanoplot 
----------------
 
Collect QC metrics on unaligned bam or fastq files with nanoplot.

pbmm2_align 
-----------
 
Uses pbmm2 to align fastq files or the raw hifi bam to the reference.

metrics_mosdepth 
----------------
 
Calculate depth stats with [Mosdepth](https://github.com/brentp/mosdepth)

set_deepvariant_regions 
-----------------------
 
Create an interval list with ScatterIntervalsByNs from GATK: [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360041416072-ScatterIntervalsByNs-Picard).
Used for creating a broken-up interval list that can be used for scattering a variant-calling pipeline in a way that will not cause problems at the edges of the intervals. 
By using large enough N blocks (so that the tools will not be able to anchor on both sides) we can be assured that the results of scattering and gathering 
the variants with the resulting interval list will be the same as calling with one large region.

deepvariant 
-----------
 
Germline variant calling with DeepVariant.

merge_filter_deepvariant 
------------------------
 
Merge deepvariant outputs, if applicable, and filter vcf.

hificnv 
-------
 
Call copy number variation and visualise results with HiFiCNV

trgt_genotyping 
---------------
 
Call tandem repeats for pathogenic and full repeats with TRGT.

sawfish 
-------
 
Call structural variants from mapped HiFi sequencing reads with Sawfish.

annotSV 
-------
 
Annotate and rank structural variants with AnnotSV.

hiphase 
-------
 
Phase variant calls with HiPhase.

report_cpsr 
-----------
 
Creates a cpsr germline report (https://sigven.github.io/cpsr/)
input: annotated/filter vcf
output: html report and addtional flat files

multiqc 
-------
 
Aggregate results from bioinformatics analyses across many samples into a single report.
MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).

