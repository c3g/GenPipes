<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Methyl-Seq Pipeline](#methyl-seq-pipeline)
  - [The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage](#the-pipeline-is-designed-to-be-run-on-a-cluster-and-is-configured-using-a-configuration-file-the-pipeline-can-be-run-in-a-single-step-or-in-multiple-steps-the-pipeline-can-also-be-run-in-parallel-to-process-multiple-samples-simultaneously%0Ausage)
  - [picard_sam_to_fastq](#picard_sam_to_fastq)
  - [trimmomatic](#trimmomatic)
  - [merge_trimmomatic_stats](#merge_trimmomatic_stats)
  - [bismark_align](#bismark_align)
  - [add_bam_umi](#add_bam_umi)
  - [sambamba_merge_sam_files](#sambamba_merge_sam_files)
  - [picard_remove_duplicates](#picard_remove_duplicates)
  - [metrics](#metrics)
  - [methylation_call](#methylation_call)
  - [wiggle_tracks](#wiggle_tracks)
  - [methylation_profile](#methylation_profile)
  - [ihec_sample_metrics_report](#ihec_sample_metrics_report)
  - [bis_snp](#bis_snp)
  - [filter_snp_cpg](#filter_snp_cpg)
  - [prepare_methylkit](#prepare_methylkit)
  - [methylkit_differential_analysis](#methylkit_differential_analysis)
  - [multiqc](#multiqc)
  - [cram_output](#cram_output)
  - [gembs_prepare](#gembs_prepare)
  - [gembs_map](#gembs_map)
  - [gembs_call](#gembs_call)
  - [gembs_bcf_to_vcf](#gembs_bcf_to_vcf)
  - [gembs_format_cpg_report](#gembs_format_cpg_report)
  - [dragen_bedgraph](#dragen_bedgraph)
  - [gembs_report](#gembs_report)
  - [dragen_align](#dragen_align)
  - [sort_dragen_sam](#sort_dragen_sam)
  - [dragen_methylation_call](#dragen_methylation_call)
  - [split_dragen_methylation_report](#split_dragen_methylation_report)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

Methyl-Seq Pipeline
================

The GenPIpes Methyl-Seq pipeline now has four protocols.
 1. bismark
 2. gembs
 3. hybrid
 4. dragen

The "bismark" protocol uses Bismark to align reads to the reference genome. Picard is used to mark and remove duplicates and generate metric files. 

The "gembs" procotol uses GemBS for mapping and methylation and variant calling (http://statgen.cnag.cat/GEMBS/UserGuide/_build/html/index.html).

The "hybrid" protocl uses [Illumina Dragen Bio-IT processor](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html) and [dragen](https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/FrontPages/DRAGEN.htm) software to align the reads to the reference genome. All the other steps are common with bismark protocol. The "dragen" protocol uses Dragen to align reads to the reference genome, call methylation, mark and remove duplicates.

Although dragen provides higher rate of mapping percentage with in a very short time duration (approximately three hours compared to 30 hours from bismark), it only accessible through McGill Genome Center cluster Abacus and The jobs cannot be submitted to any of the HPCs from the [Digital Research Aliance](https://status.computecanada.ca/). Importantly, the user needs to have permission to submit jobs to Abacus. Therefore, other users may continue to use only bismark protocol since it works in all the clusters.

However, if you would like to setup and use dragen in own cluster please refer to our [GenPipes Documentation](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_wgs_methylseq.html)

A pipeline for processing and analyzing bisulfite sequencing data. The pipeline uses Bismark to align reads and extract methylation information, and Picard to remove duplicates, add read groups and index the BAM files. The pipeline also computes metrics and generates coverage tracks per sample. The pipeline currently supports the following protocols: bismark, hybrid and dragen.

The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage
-----

```
#!text
usage: genpipes methylseq [-h] [--clean] -c CONFIG [CONFIG ...]
                          [--container {wrapper, singularity} <IMAGE PATH>]
                          [-f] [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                          [--genpipes_file GENPIPES_FILE]
                          [-j {pbs,batch,daemon,slurm}] [--json-pt]
                          [-l {debug,info,warning,error,critical}]
                          [-o OUTPUT_DIR] [--sanity-check] [-s STEPS]
                          [--wrap [WRAP]] -r READSETS_FILE [-d DESIGN_FILE]
                          [-v] [-t {bismark,gembs,hybrid,dragen}]

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
  -t, --type {bismark,gembs,hybrid,dragen}
                        Type of pipeline (default bismark)

Steps:

Protocol bismark
1 picard_sam_to_fastq
2 trimmomatic
3 merge_trimmomatic_stats
4 bismark_align
5 add_bam_umi
6 sambamba_merge_sam_files
7 picard_remove_duplicates
8 metrics
9 methylation_call
10 wiggle_tracks
11 methylation_profile
12 ihec_sample_metrics_report
13 bis_snp
14 filter_snp_cpg
15 prepare_methylkit
16 methylkit_differential_analysis
17 multiqc
18 cram_output

Protocol gembs
1 picard_sam_to_fastq
2 trimmomatic
3 merge_trimmomatic_stats
4 gembs_prepare
5 gembs_map
6 picard_remove_duplicates
7 metrics
8 gembs_call
9 gembs_bcf_to_vcf
10 gembs_format_cpg_report
11 methylation_profile
12 dragen_bedgraph
13 wiggle_tracks
14 ihec_sample_metrics_report
15 gembs_report
16 filter_snp_cpg
17 prepare_methylkit
18 methylkit_differential_analysis
19 multiqc
20 cram_output

Protocol hybrid
1 picard_sam_to_fastq
2 trimmomatic
3 merge_trimmomatic_stats
4 dragen_align
5 add_bam_umi
6 sambamba_merge_sam_files
7 picard_remove_duplicates
8 metrics
9 methylation_call
10 wiggle_tracks
11 methylation_profile
12 ihec_sample_metrics_report
13 bis_snp
14 filter_snp_cpg
15 prepare_methylkit
16 methylkit_differential_analysis
17 multiqc
18 cram_output

Protocol dragen
1 picard_sam_to_fastq
2 trimmomatic
3 merge_trimmomatic_stats
4 dragen_align
5 add_bam_umi
6 sambamba_merge_sam_files
7 sort_dragen_sam
8 metrics
9 dragen_methylation_call
10 split_dragen_methylation_report
11 methylation_profile
12 dragen_bedgraph
13 wiggle_tracks
14 ihec_sample_metrics_report
15 bis_snp
16 filter_snp_cpg
17 prepare_methylkit
18 methylkit_differential_analysis
19 multiqc
20 cram_output
```

picard_sam_to_fastq 
-------------------
 
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

trimmomatic 
-----------
 
Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', parameter 'adapter_fasta'),
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

bismark_align 
-------------
 
Align reads with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf).

add_bam_umi 
-----------
 
Add read UMI tag to individual bam files using [fgbio](https://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html).

sambamba_merge_sam_files 
------------------------
 
BAM readset files are merged into one file per sample. Merge is done using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

This step takes as input files:

1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
2. Else, BAM files from the readset file

picard_remove_duplicates 
------------------------
 
Remove duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
(for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
will be removed as a duplicate in the BAM file. Removing duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
This step is used in bismark and hybrid protocols.

metrics 
-------
 
Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads,
Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
bases which have at least 50 reads). A TDF (.tdf) coverage track is also generated at this step
for easy visualization of coverage in the IGV browser.

methylation_call 
----------------
 
The script reads in a bisulfite read alignment file produced by the Bismark bisulfite mapper
and extracts the methylation information for individual cytosines.
The methylation extractor outputs result files for cytosines in CpG, CHG and CHH context.
It also outputs bedGraph, a coverage file from positional methylation data and cytosine methylation report.

wiggle_tracks 
-------------
 
Generate wiggle tracks suitable for multiple browsers, to show coverage and methylation data.
When using GRCh37 build of Human genome, to be compatible with the UCSC Genome Browser we only keep chromosomes 1-22, X, Y and MT,
and we also rename them by prefixing "chr" to the chromosome anme (e.g. "1" becomes "chr1"), and changing the mitocondrial chromosome from "MT" to "chrM", and keeping the GRCh37 coordinates.

methylation_profile 
-------------------
 
Generation of a CpG methylation profile by combining both forward and reverse strand Cs.
Also generating of all the methylatoin metrics : CpG stats, pUC19 CpG stats, lambda conversion rate, median CpG coverage, GC bias.

ihec_sample_metrics_report 
--------------------------
 
Retrieve the computed metrics which fit the IHEC standards and build a tsv report table for IHEC.
Note: The dragen protocol does not generate a metric for estimated library size. You will have to run Picard separately for this metric.

bis_snp 
-------
 
SNP calling with [BisSNP](https://people.csail.mit.edu/dnaase/bissnp2011/).

filter_snp_cpg 
--------------
 
SNP CpGs filtering.

prepare_methylkit 
-----------------
 
Prepare input file for [methylKit](https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html) differential analysis.

methylkit_differential_analysis 
-------------------------------
 
Run methylKit to get DMCs & DMRs for different design comparisons.

multiqc 
-------
 
Aggregate results from bioinformatics analyses across many samples into a single report.
MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).

cram_output 
-----------
 
Generate long term storage version of the final alignment files in CRAM format.
Using this function will add the orginal final bam file to the removable file list.

gembs_prepare 
-------------
 
Prepare metadata and config files for mapping with gemBS.

gembs_map 
---------
 
Map reads to reference genome with GemBS's gem-mapper.

gembs_call 
----------
 
Methylation calling with bs_call as part of GemBS pipeline

gembs_bcf_to_vcf 
----------------
 
Create vcf of SNPs with bedtools intersect, by intersecting gemBS .bcf with SNP DB.

gembs_format_cpg_report 
-----------------------
 
Reformat gemBS output to match bismark and dragen output, so following steps can be followed. 

dragen_bedgraph 
---------------
 
Creates bedgraph file from combined strand CpG file

gembs_report 
------------
 
GemBS report

dragen_align 
------------
 
Align reads with [Dragen](https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/FrontPages/DRAGEN.htm) both hybrid and dragen protocols use this step to align reads.
The Dragen parameters can be changed using other_options of the ini configuration.

sort_dragen_sam 
---------------
 
Coordinate sorting the bam file resulted from dragen and create an index.

dragen_methylation_call 
-----------------------
 
Call methylation with Dragen using the 2nd run of Dragen alignment.

split_dragen_methylation_report 
-------------------------------
 
Dragen methylation report contains all three methylation context.
To create combined CSV CpGs should be extracted from the dragen methylation report.

