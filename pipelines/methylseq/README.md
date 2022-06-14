[TOC]


    Methyl-Seq Pipeline
    ================

The GenPIpes Methyl-Seq pipeline now has three protocols.
 1. bismark
 2. hybrid
 3. dragen

The "bismark" protocol uses Bismark to align reads to the reference genome.
Picard is used to mark and remove duplicates and generate metric files. The
"hybrid" protocl uses [Illumina Dragen Bio-IT processor](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html)
and [dragen](https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/FrontPages/DRAGEN.htm)
software to align the reads to the reference genome. All the other steps are common
with bismark protocol. The "dragen" protocol uses Dragen to align reads to the
reference genome, call methylation, mark and remove duplicates.

Although dragen provides higher rate of mapping percentage with in a very short time
duration (approximately three hours compared to 30 hours from bismark), it only
accessible through McGill Genome Center cluster Abacus and The jobs cannot be
submitted to any of the HPCs from the
[Digital Research Aliance](https://status.computecanada.ca/).
Importantly, the user needs to have permission to submit jobs to Abacus. Therefore,
other users may continue to use only bismark protocol since it works in all the clusters.

However, if you would like to setup and use dragen in own cluster please refer our
[GenPipes Documentation](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_wgs_methylseq.html)



Usage
-----
```
#!text

usage: methylseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                    [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                    [--no-json] [--report] [--clean]
                    [-l {debug,info,warning,error,critical}] [--sanity-check]
                    [--container {wrapper, singularity} <IMAGE PATH>]
                    [--genpipes_file GENPIPES_FILE] [-d DESIGN]
                    [-t {bismark,hybrid,dragen}] [-r READSETS] [-v]

Version: 4.3.0

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
  -d DESIGN, --design DESIGN
                        design file
  -t {bismark,hybrid,dragen}, --type {bismark,hybrid,dragen}
                        MethylSeq analysis type
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
```
![methylseq bismark workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_methylseq_bismark.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_methylseq_bismark.png)
```
bismark:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- bismark_align
5- add_bam_umi
6- sambamba_merge_sam_files
7- picard_remove_duplicates
8- metrics
9- methylation_call
10- wiggle_tracks
11- methylation_profile
12- ihec_sample_metrics_report
13- bis_snp
14- filter_snp_cpg
15- prepare_methylkit
16- cram_output
----
```
![methylseq hybrid workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_methylseq_hybrid.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_methylseq_hybrid.png)
```
hybrid:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- dragen_align
5- add_bam_umi
6- sambamba_merge_sam_files
7- picard_remove_duplicates
8- metrics
9- methylation_call
10- wiggle_tracks
11- methylation_profile
12- ihec_sample_metrics_report
13- bis_snp
14- filter_snp_cpg
15- prepare_methylkit
16- cram_output
----
```
![methylseq dragen workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_methylseq_dragen.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_methylseq_dragen.png)
```
dragen:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- dragen_align
5- add_bam_umi
6- sambamba_merge_sam_files
7- sort_dragen_sam
8- metrics
9- dragen_methylation_call
10- split_dragen_methylation_report
11- methylation_profile
12- dragen_bedgraph
13- wiggle_tracks
14- ihec_sample_metrics_report
15- bis_snp
16- filter_snp_cpg
17- prepare_methylkit
18- cram_output

```

picard_sam_to_fastq
-------------------
Converts SAM/BAM files from the input readset file into FASTQ format.
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

cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format.
Using this function will include the orginal final bam file into the  removable file list.

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


