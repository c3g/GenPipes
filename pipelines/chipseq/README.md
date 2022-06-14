[TOC]


ChIP-Seq Pipeline
=================

ChIP-Seq experiments allows the Isolation and sequencing of genomic DNA bound by a specific transcription factor,
covalently modified histone, or other nuclear protein. The pipeline starts by trimming adaptors and
low quality bases and mapping the reads (single end or paired end ) to a reference genome using bwa.
Reads are filtered by mapping quality and duplicate reads are marked. Then, Homer quality control routines
are used to provide information and feedback about the quality of the experiment. Peak calls is executed by MACS
and annotation and motif discovery for narrow peaks are executed using Homer. Statistics of annotated peaks
are produced for narrow peaks and a standard report is generated.

An example of the ChIP-Seq report for an analysis on public ENCODE data is available for illustration purpose only:
[ChIP-Seq report](http://gqinnovationcenter.com/services/bioinformatics/tools/chipReport/index.html).

[Here](https://bitbucket.org/mugqic/mugqic_pipelines/downloads/MUGQIC_Bioinfo_ChIP-Seq.pptx)
is more information about ChIP-Seq pipeline that you may find interesting.


Usage
-----
```
#!text

usage: chipseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                  [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                  [--no-json] [--report] [--clean]
                  [-l {debug,info,warning,error,critical}] [--sanity-check]
                  [--container {wrapper, singularity} <IMAGE PATH>]
                  [--genpipes_file GENPIPES_FILE] [-d DESIGN]
                  [-t {chipseq,atacseq}] [-r READSETS] [-v]

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
  -t {chipseq,atacseq}, --type {chipseq,atacseq}
                        Type of pipeline (default chipseq)
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
```
![chipseq chipseq workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_chipseq_chipseq.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_chipseq_chipseq.png)
```
chipseq:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- mapping_bwa_mem_sambamba
5- sambamba_merge_bam_files
6- sambamba_mark_duplicates
7- sambamba_view_filter
8- metrics
9- homer_make_tag_directory
10- qc_metrics
11- homer_make_ucsc_file
12- macs2_callpeak
13- homer_annotate_peaks
14- homer_find_motifs_genome
15- annotation_graphs
16- run_spp
17- differential_binding
18- ihec_metrics
19- multiqc_report
20- cram_output
21- gatk_haplotype_caller
22- merge_and_call_individual_gvcf
----
```
![chipseq atacseq workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_chipseq_atacseq.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_chipseq_atacseq.png)
```
atacseq:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- mapping_bwa_mem_sambamba
5- sambamba_merge_bam_files
6- sambamba_mark_duplicates
7- sambamba_view_filter
8- metrics
9- homer_make_tag_directory
10- qc_metrics
11- homer_make_ucsc_file
12- macs2_atacseq_callpeak
13- homer_annotate_peaks
14- homer_find_motifs_genome
15- annotation_graphs
16- run_spp
17- differential_binding
18- ihec_metrics
19- multiqc_report
20- cram_output
21- gatk_haplotype_caller
22- merge_and_call_individual_gvcf

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

mapping_bwa_mem_sambamba
------------------------
The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

This step takes as input files:

1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

sambamba_merge_bam_files
------------------------
BAM readset files are merged into one file per sample. Merge is done using [Sambamba]().

This step takes as input files:

1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
2. Else, BAM files from the readset file

sambamba_mark_duplicates
------------------------
Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
(for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
will be marked as a duplicate in the BAM file. Marking duplicates is done using [Sambamba]().

sambamba_view_filter
--------------------
Filter out unmapped reads and low quality reads [Sambamba](http://www.htslib.org/).

metrics
-------
The number of raw/filtered and aligned reads per sample are computed at this stage.

homer_make_tag_directory
------------------------
The Homer Tag directories, used to check for quality metrics, are computed at this step.

qc_metrics
----------
Sequencing quality metrics as tag count, tag autocorrelation, sequence bias and GC bias are generated.

homer_make_ucsc_file
--------------------
Wiggle Track Format files are generated from the aligned reads using Homer.
The resulting files can be loaded in browsers like IGV or UCSC.

macs2_callpeak
--------------
Peaks are called using the MACS2 software. Different calling strategies are used for narrow and broad peaks.
The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run.
The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100.
The default mfold parameter of MACS2 is [10,30].

homer_annotate_peaks
--------------------
The peaks called previously are annotated with HOMER using RefSeq annotations for the reference genome.
Gene ontology and genome ontology analysis are also performed at this stage.

homer_find_motifs_genome
------------------------
De novo and known motif analysis per design are performed using HOMER.

annotation_graphs
-----------------
The peak location statistics. The following peak location statistics are generated per design:
proportions of the genomic locations of the peaks. The locations are: Gene (exon or intron),
Proximal ([0;2] kb upstream of a transcription start site), Distal ([2;10] kb upstream
of a transcription start site), 5d ([10;100] kb upstream of a transcription start site),
Gene desert (>= 100 kb upstream or downstream of a transcription start site), Other (anything
not included in the above categories); The distribution of peaks found within exons and introns;
The distribution of peak distance relative to the transcription start sites (TSS);
the Location of peaks per design.

run_spp
-------
runs spp to estimate NSC and RSC ENCODE metrics. For more information: https://github.com/kundajelab/phantompeakqualtools

differential_binding
--------------------
Performs differential binding analysis using [DiffBind](http://bioconductor.org/packages/release/bioc/html/DESeq.html)
Merge the results of the analysis in a single csv file.
html report will be generated to QC samples and check how well differential binding analysis was performed.

ihec_metrics
------------
Generate IHEC's standard metrics.

multiqc_report
--------------
A quality control report for all samples is generated.
For more detailed information about the MultiQc visit: [MultiQc] (http://multiqc.info/)

cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format
Using this function will include the orginal final bam file into the  removable file list

gatk_haplotype_caller
---------------------
GATK haplotype caller for snps and small indels.

merge_and_call_individual_gvcf
------------------------------
Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.

macs2_atacseq_callpeak
----------------------
Peaks are called using the MACS2 software. Different calling strategies are used for narrow and broad peaks.
The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run.
The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100.
The default mfold parameter of MACS2 is [10,30].


