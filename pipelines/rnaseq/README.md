[TOC]



Usage
-----
```
#!text

usage: rnaseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                 [--no-json] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [--sanity-check]
                 [--container {wrapper, singularity} <IMAGE PATH>]
                 [--genpipes_file GENPIPES_FILE] [-t {stringtie,cufflinks}]
                 [-d DESIGN] [-r READSETS] [-v]

Version: 4.1.2

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

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
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -t {stringtie,cufflinks}, --type {stringtie,cufflinks}
                        RNAseq analysis type
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
```
![rnaseq stringtie workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq_stringtie.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq_stringtie.png)
```
stringtie:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- star
5- picard_merge_sam_files
6- picard_sort_sam
7- picard_mark_duplicates
8- picard_rna_metrics
9- estimate_ribosomal_rna
10- bam_hard_clip
11- rnaseqc
12- wiggle
13- raw_counts
14- raw_counts_metrics
15- stringtie
16- stringtie_merge
17- stringtie_abund
18- ballgown
19- differential_expression
20- cram_output
----
```
![rnaseq cufflinks workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq_cufflinks.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq_cufflinks.png)
```
cufflinks:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- star
5- picard_merge_sam_files
6- picard_sort_sam
7- picard_mark_duplicates
8- picard_rna_metrics
9- estimate_ribosomal_rna
10- bam_hard_clip
11- rnaseqc
12- wiggle
13- raw_counts
14- raw_counts_metrics
15- cufflinks
16- cuffmerge
17- cuffquant
18- cuffdiff
19- cuffnorm
20- fpkm_correlation_matrix
21- gq_seq_utils_exploratory_analysis_rnaseq
22- differential_expression
23- differential_expression_goseq
24- ihec_metrics
25- cram_output

```
picard_sam_to_fastq
-------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

trimmomatic
-----------
Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
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

star
----
The filtered reads are aligned to a reference genome. The alignment is done per readset of sequencing
using the [STAR](https://code.google.com/p/rna-star/) software. It generates a Binary Alignment Map file (.bam).

This step takes as input files:

1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

picard_merge_sam_files
----------------------
BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

picard_sort_sam
---------------
The alignment file is reordered (QueryName) using [Picard](http://broadinstitute.github.io/picard/). The QueryName-sorted bam files will be used to determine raw read counts.

picard_mark_duplicates
----------------------
Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
(for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).

picard_rna_metrics
------------------
Computes a series of quality control metrics using both CollectRnaSeqMetrics and CollectAlignmentSummaryMetrics functions
metrics are collected using [Picard](http://broadinstitute.github.io/picard/).

estimate_ribosomal_rna
----------------------
Use bwa mem to align reads on the rRNA reference fasta and count the number of read mapped
The filtered reads are aligned to a reference fasta file of ribosomal sequence. The alignment is done per sequencing readset.
The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
BWA output BAM files are then sorted by coordinate using [Picard](http://broadinstitute.github.io/picard/).

This step takes as input files:

readset Bam files

bam_hard_clip
-------------
Generate a hardclipped version of the bam for the toxedo suite which doesn't support this official sam feature.

rnaseqc
-------
Computes a series of quality control metrics using [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc).

wiggle
------
Generate wiggle tracks suitable for multiple browsers.

raw_counts
----------
Count reads in features using [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html).

raw_counts_metrics
------------------
Create rawcount matrix, zip the wiggle tracks and create the saturation plots based on standardized read counts.

stringtie
---------
Assemble transcriptome using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
Warning: Still in testing.

stringtie_merge
---------------
Merge assemblies into a master teranscriptome reference using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
Warning: still in testing

stringtie_abund
---------------
Assemble transcriptome and compute RNA-seq expression using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
Warning: Still in testing.

ballgown
--------
[Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is used to calculate differential transcript and gene expression levels and test them for significant differences.

Warning: still in testing

differential_expression
-----------------------
Performs differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
Merge the results of the analysis in a single csv file.

cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format
Using this function will include the orginal final bam file into the  removable file list

cufflinks
---------
Compute RNA-Seq data expression using [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/).
Warning: It needs to use a hard clipped bam file while Tuxedo tools do not support official soft clip SAM format

cuffmerge
---------
Merge assemblies into a master transcriptome reference using [cuffmerge](http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/).

cuffquant
---------
Compute expression profiles (abundances.cxb) using [cuffquant](http://cole-trapnell-lab.github.io/cufflinks/cuffquant/).
Warning: It needs to use a hard clipped bam file while Tuxedo tools do not support official soft clip SAM format

cuffdiff
--------
[Cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/) is used to calculate differential transcript expression levels and test them for significant differences.

cuffnorm
--------
Global normalization of RNA-Seq expression levels using [Cuffnorm](http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/).

fpkm_correlation_matrix
-----------------------
Compute the pearson corrleation matrix of gene and transcripts FPKM. FPKM data are those estimated by cuffnorm.

gq_seq_utils_exploratory_analysis_rnaseq
----------------------------------------
Exploratory analysis using the gqSeqUtils R package.

differential_expression_goseq
-----------------------------
Gene Ontology analysis for RNA-Seq using the Bioconductor's R package [goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html).
Generates GO annotations for differential gene expression analysis.

ihec_metrics
------------
Generate IHEC's standard metrics.


