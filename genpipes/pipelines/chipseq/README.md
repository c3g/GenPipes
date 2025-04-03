<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [ChIP-Seq Pipeline](#chip-seq-pipeline)
  - [The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage](#the-pipeline-is-designed-to-be-run-on-a-cluster-and-is-configured-using-a-configuration-file-the-pipeline-can-be-run-in-a-single-step-or-in-multiple-steps-the-pipeline-can-also-be-run-in-parallel-to-process-multiple-samples-simultaneously%0Ausage)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

ChIP-Seq Pipeline
=================

A pipeline to process ChIP-seq data. The pipeline is designed to handle both paired-end and single-end reads and can be used to process data from any Illumina sequencer. The pipeline uses Trimmomatic to trim reads and remove Illumina adapters, BWA to align reads to the reference genome, Sambamba to sort and filter BAM files, and Picard to mark duplicates and collect quality metrics. The pipeline also uses MACS2 to call peaks, HOMER to annotate peaks, and DiffBind to perform differential binding analysis.

The pipeline takes as input a readset file and a design file. The readset file contains the list of samples and readsets, while the design file contains the list of contrasts to be analyzed. The pipeline outputs BAM files, peak calls, and differential binding results. The pipeline also generates quality metrics and reports for each sample.

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

