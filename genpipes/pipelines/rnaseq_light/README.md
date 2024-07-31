[TOC]

Rnaseq_light Pipeline
================

Usage
-----


```
#!text

usage: genpipes rnaseq_light [-h] -c CONFIG [CONFIG ...] [-s STEPS]
                             [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}]
                             [-f] [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                             [--no-json] [--json-pt] [--clean]
                             [--container {wrapper, singularity} <IMAGE PATH>]
                             [--genpipes_file GENPIPES_FILE]
                             [-l {debug,info,warning,error,critical}]
                             [--sanity-check] [--wrap [WRAP]] -r READSETS_FILE
                             [-d DESIGN_FILE] [-v]

Version: 5.0.0

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

options:
  -h, --help            show this help message and exit
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
  --force_mem_per_cpu FORCE_MEM_PER_CPU
                        Take the mem input in the ini file and force to have a
                        minimum of mem_per_cpu by correcting the number of cpu
                        (default: None)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a validsingularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  --wrap [WRAP]         Path to the genpipe cvmfs wrapper script. Default is g
                        enpipes/ressources/container/bin/container_wrapper.sh.
                        This is a convenience options for using genpipes in a
                        container
  -r READSETS_FILE, --readsets READSETS_FILE
                        readset file
  -d DESIGN_FILE, --design DESIGN_FILE
                        design file
  -v, --version         show the version information and exit

Summary:

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
    in the differential analyses. The design file format is described
    [here](https://genpipes.readthedocs.io/en/latest/get-started/concepts/design_file.html). In addition,
    [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is used to calculate differential
    transcript and gene expression levels and test them for significant differences.
    It can also take a batch file (optional) which will be used to correct for batch effects
    in the differential analyses. The batch file format is described
    [here](https://bitbucket.org/mugqic/mugqic_pipelines/src#markdown-header-batch-file)

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

    
Steps:

Protocol default
0 picard_sam_to_fastq
1 trimmomatic
2 merge_trimmomatic_stats
3 kallisto
4 kallisto_count_matrix
5 gq_seq_utils_exploratory_analysis_rnaseq_light
6 sleuth_differential_expression
7 multiqc
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

kallisto 
--------
 
Run [Kallisto](https://pachterlab.github.io/kallisto/about.html) on fastq files for a fast esimate of abundance.

kallisto_count_matrix 
---------------------
 
Use the output from Kallisto to create a transcript count matrix.
Create a summary table to be included in the multiqc report.

gq_seq_utils_exploratory_analysis_rnaseq_light 
----------------------------------------------
 
Exploratory analysis using the gqSeqUtils R package adapted for RnaSeqLight.

sleuth_differential_expression 
------------------------------
 
Performs differential gene expression analysis using [Sleuth](http://pachterlab.github.io/sleuth/).
Analysis are performed both at a transcript and gene level, using two different tests: LRT and WT.

multiqc 
-------
 
A quality control report for all samples is generated.
For more detailed information about MultiQC visit: [MultiQC](http://multiqc.info/)

