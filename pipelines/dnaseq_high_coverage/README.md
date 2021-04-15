[TOC]


DnaSeq high Coverage Pipeline
=================

The DnaSeq high Coverage Pipeline is based on the DNA-Seq pipeline and follow the first initial steps.
The difference starts at the Mark Duplicates step. Since this data is high coverage Mark Dup is not run.
Recalibration is not run either because typically, these datasets are targetted with amplicons or custom
capture which render recalibration useless.

Also variant calling is done only using frequency. Not Bayesian callers are used because these typically
don't fare well with the high coverage.


Usage
-----
```
#!text

usage: dnaseq_high_coverage.py [-h] [--help] [-c CONFIG [CONFIG ...]]
                               [-s STEPS] [-o OUTPUT_DIR]
                               [-j {pbs,batch,daemon,slurm}] [-f] [--no-json]
                               [--report] [--clean]
                               [-l {debug,info,warning,error,critical}]
                               [--sanity-check]
                               [--container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}]
                               [-t {mugqic,mpileup,light}] [-r READSETS] [-v]

Version: 3.1.5

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

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
  --container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}
                        run pipeline inside a container providing a container
                        image path or accessible docker/singularity hub path
  -t {mugqic,mpileup,light}, --type {mugqic,mpileup,light}
                        DNAseq analysis type
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
```
![dnaseq_high_coverage workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_high_coverage.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_high_coverage.png)
```
------
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- bwa_mem_picard_sort_sam
5- picard_merge_sam_files
6- gatk_indel_realigner
7- merge_realigned
8- picard_fixmate
9- metrics
10- picard_calculate_hs_metrics
11- gatk_callable_loci
12- call_variants
13- preprocess_vcf
14- snp_effect
15- gemini_annotations
16- cram_output

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

bwa_mem_picard_sort_sam
-----------------------
The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
BWA output BAM files are then sorted by coordinate using [Picard](http://broadinstitute.github.io/picard/).

This step takes as input files:

1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

picard_merge_sam_files
----------------------
BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

This step takes as input files:

1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
2. Else, BAM files from the readset file

gatk_indel_realigner
--------------------
Insertion and deletion realignment is performed on regions where multiple base mismatches
are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
Such regions will introduce false positive variant calls which may be filtered out by realigning
those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
The reference genome is divided by a number regions given by the `nb_jobs` parameter.

merge_realigned
---------------
BAM files of regions of realigned reads are merged per sample using [Picard](http://broadinstitute.github.io/picard/).

picard_fixmate
--------------
Verify mate-pair information between mates and fix if needed.
This ensures that all mate-pair information is in sync between each read and its mate pair.
Fix is done using [Picard](http://broadinstitute.github.io/picard/).

metrics
-------
Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads,
Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
bases which have at least 50 reads). A TDF (.tdf) coverage track is also generated at this step
for easy visualization of coverage in the IGV browser.

picard_calculate_hs_metrics
---------------------------
Compute on target percent of hybridisation based capture.

gatk_callable_loci
------------------
Computes the callable region or the genome as a bed track.

call_variants
-------------
VarScan caller for insertions and deletions.

preprocess_vcf
--------------
Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html
Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and 
vcf FORMAT modification for correct loading into gemini

snp_effect
----------
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

gemini_annotations
------------------
Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html

cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format
Using this function will include the orginal final bam file into the  removable file list 


