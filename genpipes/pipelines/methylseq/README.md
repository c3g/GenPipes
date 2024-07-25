[TOC]

Methylseq Pipeline
================

Usage
-----


```
#!text

usage: genpipes methylseq [-h] -c CONFIG [CONFIG ...] [-s STEPS]
                          [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                          [--force_mem_per_cpu FORCE_MEM_PER_CPU] [--no-json]
                          [--json-pt] [--clean]
                          [--container {wrapper, singularity} <IMAGE PATH>]
                          [--genpipes_file GENPIPES_FILE]
                          [-l {debug,info,warning,error,critical}]
                          [--sanity-check] [--wrap [WRAP]] -r READSETS_FILE
                          [-d DESIGN_FILE] [-v] [-t {bismark,hybrid,dragen}]

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
  -t {bismark,hybrid,dragen}, --type {bismark,hybrid,dragen}
                        Type of pipeline (default chipseq)

Steps:

Protocol bismark
0 picard_sam_to_fastq
1 trimmomatic
2 merge_trimmomatic_stats
3 bismark_align
4 add_bam_umi
5 sambamba_merge_sam_files
6 picard_remove_duplicates
7 metrics
8 methylation_call
9 wiggle_tracks
10 methylation_profile
11 ihec_sample_metrics_report
12 bis_snp
13 filter_snp_cpg
14 prepare_methylkit
15 methylkit_differential_analysis
16 multiqc
17 cram_output

Protocol hybrid
0 picard_sam_to_fastq
1 trimmomatic
2 merge_trimmomatic_stats
3 dragen_align
4 add_bam_umi
5 sambamba_merge_sam_files
6 picard_remove_duplicates
7 metrics
8 methylation_call
9 wiggle_tracks
10 methylation_profile
11 ihec_sample_metrics_report
12 bis_snp
13 filter_snp_cpg
14 prepare_methylkit
15 methylkit_differential_analysis
16 multiqc
17 cram_output

Protocol dragen
0 picard_sam_to_fastq
1 trimmomatic
2 merge_trimmomatic_stats
3 dragen_align
4 add_bam_umi
5 sambamba_merge_sam_files
6 sort_dragen_sam
7 metrics
8 dragen_methylation_call
9 split_dragen_methylation_report
10 methylation_profile
11 dragen_bedgraph
12 wiggle_tracks
13 ihec_sample_metrics_report
14 bis_snp
15 filter_snp_cpg
16 prepare_methylkit
17 methylkit_differential_analysis
18 multiqc
19 cram_output
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

dragen_bedgraph 
---------------
 
Creates bedgraph file from combined strand CpG file

