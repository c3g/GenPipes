<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [RNA-Seq Light Pipeline](#rna-seq-light-pipeline)
  - [Usage](#usage)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

RNA-Seq Light Pipeline
================

The RNA-Seq Light Pipeline is a lightweight analysis of gene expression in RNA sequencing data. 
The pipeline is based on [Kallisto](https://pachterlab.github.io/kallisto/about.html) and differential expression analysis is performed by [Sleuth](http://pachterlab.github.io/sleuth/).
It is especially useful for quick Quality Control (QC) in gene sequencing studies.
    
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

