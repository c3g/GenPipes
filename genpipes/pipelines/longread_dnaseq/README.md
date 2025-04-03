<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [LongRead DNA-Seq Pipeline](#longread-dna-seq-pipeline)
  - [Usage](#usage)

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

