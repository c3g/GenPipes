<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Methyl-Seq Pipeline](#methyl-seq-pipeline)
  - [The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage](#the-pipeline-is-designed-to-be-run-on-a-cluster-and-is-configured-using-a-configuration-file-the-pipeline-can-be-run-in-a-single-step-or-in-multiple-steps-the-pipeline-can-also-be-run-in-parallel-to-process-multiple-samples-simultaneously%0Ausage)

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
usage: genpipes [-h] [-v] [-s {bash,zsh,tcsh}]
                {ampliconseq,chipseq,covseq,dnaseq,longread_dnaseq,methylseq,nanopore_covseq,rnaseq,rnaseq_denovo_assembly,rnaseq_light,tools} ...

GenPipes consists of Python scripts which create a list of jobs running Bash
commands. Those scripts support dependencies between jobs and smart restart
mechanism if some jobs fail during pipeline execution. Jobs can be submitted
in different ways: by being sent to a PBS or a SLURM scheduler or by being run
as a series of commands in batch through a Bash script. Job commands and
parameters can be modified through several configuration files.

positional arguments:
  {ampliconseq,chipseq,covseq,dnaseq,longread_dnaseq,methylseq,nanopore_covseq,rnaseq,rnaseq_denovo_assembly,rnaseq_light,tools}
    ampliconseq         AmpliconSeq pipeline
    chipseq             ChipSeq pipeline
    covseq              CoVSeq pipeline
    dnaseq              DnaSeq pipeline
    longread_dnaseq     LongRead DnaSeq pipeline
    methylseq           MethylSeq pipeline
    nanopore_covseq     Nanopore CoVSeq pipeline
    rnaseq              RnaSeq pipeline
    rnaseq_denovo_assembly
                        RnaSeq DeNovo Assembly pipeline
    rnaseq_light        RnaSeq Light pipeline
    tools               GenPipes companion tools

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -s, --print-completion {bash,zsh,tcsh}
                        print shell completion script

