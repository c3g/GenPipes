[TOC]


Capture Hi-C Pipeline
=====================

    Hi-C experiments allow researchers to understand chromosomal folding and structure using proximity ligation techniques.
    Capture Hi-C experiments allow researchers to focus on regions of interest by capturing those regions on a chip.
    The pipeline starts by trimming adaptors and low quality bases. It then maps the reads to a reference genome using HiCUP.
    HiCUP first truncates chimeric reads, maps reads to the genome, then it filters Hi-C artifacts and removes read duplicates.
    Samples from different lanes are merged. CHiCAGO is then used to filter capture-specific artifacts and call significant 
    interactions. This pipeline also identifies enrichement of regulatory features when provided with ChIP-Seq marks and also 
    identifies enrichement of GWAS SNPs.


Usage
-----
```
#!text

usage: hicseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS] [-e Restriction_Enzyme {DpnII, MboI, HindIII, NcoI}]
                  [-o OUTPUT_DIR] [-j {pbs,batch}] [-f] [--report] [--clean]
                  [-l {debug,info,warning,error,critical}] [-d DESIGN]
                  [-r READSETS] [-v]

Version: 1.0.0

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -e Restriction_Enzyme --enzyme Restriction_Enzyme
                        restriction enzyme used in generating Hi-C library
                        enzymes currently supported: DpnII, MboI, HindIII, NcoI
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  -j {pbs,batch}, --job-scheduler {pbs,batch}
                        job scheduler type (default: pbs)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
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
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
1- samtools_bam_sort
2- picard_sam_to_fastq
3- trimmomatic
4- merge_trimmomatic_stats
5- fastq_readName_Edit
6- hicup_align
7- samtools_merge_bams

15- multiqc_report

```
1- samtools_bam_sort
----------------------
Sorts bam by readname prior to picard_sam_to_fastq step in order to minimize memory consumption.
If bam file is small and the memory requirements are reasonable, this step can be skipped.


2- picard_sam_to_fastq
----------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

3- trimmomatic
--------------
Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
only Adapter1 is used and left unchanged.

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

4- merge_trimmomatic_stats
--------------------------
The trim statistics per readset are merged at this step.


5- fastq_readName_Edit
----------------------
Removes the added /1 and /2 by picard's sam_to_fastq transformation to avoid issues with downstream software like HOMER
      

6- hicup_align
--------------------------
Paired-end Hi-C reads are truncated, mapped and filtered using HiCUP. The resulting bam file is filtered for Hi-C artifacts and duplicated reads. It is ready for use as input for downstream analysis.

For more detailed information about the HICUP process visit: [HiCUP] (https://www.bioinformatics.babraham.ac.uk/projects/hicup/overview/)
   

7- samtools_merge_bams
-----------------------
BAM readset files are merged into one file per sample. Merge is done using [samtools](http://samtools.sourceforge.net/). This step takes as input files the aligned bams/sams from the hicup_align step.
    


15- multiqc_report
-------------------
A quality control report for all samples is generated.

For more detailed information about the MultiQc visit: [MultiQc] (http://multiqc.info/)

