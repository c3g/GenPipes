<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [RNA-Seq Pipeline](#rna-seq-pipeline)
  - [Usage](#usage)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

RNA-Seq Pipeline
================

The standard MUGQIC RNA-Seq pipeline has three protocols (stringtie, variants, cancer), stringtie is the
default protocol and applicable in most cases.

All three protocols are based on the use of the [STAR aligner](https://code.google.com/p/rna-star/)
to align reads to the reference genome. These alignments are used during
downstream analysis (for example in stringtie protocol, to determine differential expression of genes and transcripts).

The [StringTie](https://ccb.jhu.edu/software/stringtie/) suite is used for differential transcript expression (DTE) analysis, whereas
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and
[edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) are used for the differential gene expression (DGE) analysis.

The "stringtie" protocol requires a design file which will be used to define the comparison groups
in the differential analyses. The design file format is described [here](https://genpipes.readthedocs.io/en/latest/get-started/concepts/design_file.html).
In addition, [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is used to calculate differential
transcript and gene expression levels and test them for significant differences.
It can also take a batch file (optional) which will be used to correct for batch effects
in the differential analyses. The batch file format is described [here](https://github.com/c3g/GenPipes?tab=readme-ov-file#batch-file)

The variants protocol is used when variant detection, is the main goal of the analysis. GATK best practices workflow
is used for variant calling. It also contains a step for annotating genes using [gemini](https://gemini.readthedocs.io/en/latest/)

The cancer protocol contains all the steps in the variants protocol but it is
specifically designed for cancer data sets due to the
complexity of cancer samples and additional analyses those projects often entail.
The goal of the cancer protocol is comparing expression to known benchmarks. In addition to the steps in the variants
protocol, it contains four specific steps. Three of them (run_star_fusion, run_arriba, run_annofuse)
are related to detection and annotation of gene fusions. For that, [Star-fusion](https://github.com/STAR-Fusion/STAR-Fusion-Tutorial/wiki),
[Arriba](https://arriba.readthedocs.io/en/latest/) and [anno-Fuse](https://rdrr.io/github/d3b-center/annoFuse/) are
used. Furthermore, [RSeQC](http://rseqc.sourceforge.net/) provides RNA-seq quality control metrics to asses the
quality of the data.
    
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

