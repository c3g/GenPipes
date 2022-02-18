[TOC]


MGI Run Processing Pipeline
================================

The standard Run Processing pipeline handles both Illumina and MGI sequencing technologies.
It uses the Illumina bcl2fastq software to convert and demultiplex Illumina base call files to fastq files.
In the case of MGI run processing, it uses fastq files produced by the MGI-G400 sequencer, or MGI-T7 base call files,
then does demultiplexing. Finally, the pipeline runs some QCs on the raw data, on the fastq and on the alignment.

Sample Sheets
-------------

The pipeline uses one input sample sheet, a tsv file having the following columns:

- ProjectLUID
- ProjectName
- ContainerName
- Position
- Index
- LibraryLUID
- LibraryProcess
- SampleLUID
- SampleName
- Reference
- Start Date
- Sample Tag
- Target Cells
- Species
- UDF/Genome Size (Mb)
- Gender
- Pool Fraction
- Capture Type
- Capture Name
- Capture REF_BED
- Library Size
- Library Kit Name
- Capture Kit Type
- Capture Bait Version
- ChIP-Seq Mark

Example:
ProjectLUID	ProjectName	ContainerName	Position	Index	LibraryLUID	LibraryProcess	SampleName	Reference	Start Date	Sample Tag	Target Cells	Species	UDF/Genome Size (Mb)	Gender	Pool Fraction	Capture Type	Capture Name	Capture REF_BED	Library Size	Library Kit Name	Capture Kit Type	Capture Bait Version	ChIP-Seq Mark
AUL208	NA_Control	V300096783	4:1	MGI09_A02_Barcode_41	2-1981727	RNASeq MGI	RNA_GM12878_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
AUL208	NA_Control	V300096783	4:1	MGI13_E02_Barcode_45	2-1981728	RNASeq MGI	RNA_GM12878_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
AUL208	NA_Control	V300096783	2:1	MGI10_B02_Barcode_42	2-1981725	RNASeq MGI	RNA_Mother_HG004_GM24143_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
AUL208	NA_Control	V300096783	2:1	MGI14_F02_Barcode_46	2-1981726	RNASeq MGI	RNA_Mother_HG004_GM24143_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
AUL208	NA_Control	V300096783	4:1	MGI10_B02_Barcode_42	2-1981725	RNASeq MGI	RNA_Mother_HG004_GM24143_1	N/A	2021-12-09	N/A	N/A	Eukaryota:Homo sapiens (Taxon ID:9606)	3257.32	F	0.125	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A


Usage
-----
```
#!text

usage: run_processing.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                         [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}]
                         [--force] [--no-json] [--report] [--clean]
                         [-l {debug,info,warning,error,critical}]
                         [--sanity-check]
                         [--container {wrapper, singularity} <IMAGE PATH>]
                         [--genpipes_file GENPIPES_FILE]
                         [-t {illumina,mgig400,mgit7}] [-r READSETS]
                         [-d RUN_DIR] [--run-id RUN_ID] [-f RAW_FLAG_DIR]
                         [--splitbarcode-demux] [--lane LANE_NUMBER]
                         [-x FIRST_INDEX] [-y LAST_INDEX]
                         [-m NUMBER_OF_MISMATCHES] [--allow-barcode-collision]
                         [-v]

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
  --force               force creation of jobs even if up to date (default:
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
  -t {illumina,mgig400,mgit7}, --type {illumina,mgig400,mgit7}
                        Sequencing technology : Illumina, MGI G400 or MGI T7
                        (mandatory)
  -r READSETS, --readsets READSETS
                        Sample sheet for the MGI run to process (mandatory)
  -d RUN_DIR, --run RUN_DIR
                        Run directory (mandatory)
  --run-id RUN_ID       Run ID. Default is parsed from the run folder
  -f RAW_FLAG_DIR, --flag RAW_FLAG_DIR
                        T7 flag files directory (mandatory for MGI T7 runs)
  --splitbarcode-demux  demultiplexing done while basecalling with MGI
                        splitBarcode (only affect MGI G400 or T7 runs)
  --lane LANE_NUMBER    Lane number (to only process the given lane)
  -x FIRST_INDEX        First index base to use for demultiplexing
                        (inclusive). The index from the sample sheet will be
                        adjusted according to that value.
  -y LAST_INDEX         Last index base to use for demultiplexing (inclusive)
  -m NUMBER_OF_MISMATCHES
                        Number of index mistmaches allowed for demultiplexing
                        (default 1). Barcode collisions are always checked.
  --allow-barcode-collision
                        Allow barcode collision by not comparing barcode
                        sequences to each other (usually decreases the
                        demultiplexing efficiency).
  -v, --version         show the version information and exit

Steps:
------

----
illumina:
1- index
2- fastq
3- qc_graphs
4- fastp
5- fastqc
6- blast
7- align
8- picard_mark_duplicates
9- metrics
10- md5
11- report
12- copy
13- final_notification
----
mgig400:
1- fastq
2- qc_graphs
3- fastp
4- fastqc
5- blast
6- align
7- picard_mark_duplicates
8- metrics
9- md5
10- report
11- copy
12- final_notification
----
mgit7:
1- basecall
2- fastq
3- qc_graphs
4- fastp
5- fastqc
6- blast
7- align
8- picard_mark_duplicates
9- metrics
10- md5
11- report
12- copy
13- final_notification

```
index
-----
Generate a file with all the indexes found in the index-reads of the run.

The input barcode file is a two columns tsv file. Each line has a
`barcode_sequence` and the corresponding `barcode_name`. This file can be
generated by a LIMS.

The output is a tsv file named `RUNFOLDER_LANENUMBER.metrics` that will be
saved in the output directory. This file has four columns, the barcode/index
sequence, the index name, the number of reads and the number of reads that have
passed the filter.

fastq
-----
For Illumina

    Launch fastq generation from Illumina raw data using BCL2FASTQ conversion
    software.

    The index base mask is calculated according to the sample and run configuration;
    and also according the mask parameters received (first/last index bases). The
    Casava sample sheet is generated with this mask. The default number of
    mismatches allowed in the index sequence is 1 and can be overrided with an
    command line argument. Demultiplexing always occurs even when there is only one
    sample in the lane, because then we merge undertermined reads.

    An optional notification command can be launched to notify the start of the
    fastq generation with the calculated mask.

For MGI-G400

    First copy all the files of the lane from the sequencer deposit folder
    to the processing folder, into "raw_fastq".
    Then, perform demultplexing of the reads with fgbio DemuxFastqs

For MGI-T7

    Perform demultiplexing of the reads with fgbio DemuxFastqs
    (skipped with --splitbarcode-demux)

qc_graphs
---------
Generate some QC Graphics and a summary XML file for each sample using
[BVATools](https://bitbucket.org/mugqic/bvatools/).

Files are created in a 'qc' subfolder of the fastq directory. Examples of
output graphic:

- Per cycle qualities, sequence content and sequence length;
- Known sequences (adaptors);
- Abundant Duplicates;

fastp
-----
Generate basic QC metrics.

fastqc
------

blast
-----
Run blast on a subsample of the reads of each sample to find the 20 most
frequent hits.

The `runBlast.sh` tool from MUGQIC Tools is used. The number of reads to
subsample can be configured by sample or for the whole lane. The output will be
in the `Blast_sample` folder, under the Unaligned folder.

align
-----
Align the reads from the fastq file, sort the resulting .bam and create an index
of that .bam.

An basic aligment is performed on a sample when the `SampleRef` field of the
MGI sample sheet match one of the regexp in the configuration file and the
corresponding genome (and indexes) are installed.

`STAR` is used as a splice-junctions aware aligner when the sample
`library_source` is `cDNA` or contains `RNA`; otherwise `BWA_mem` is used to
align the reads.

picard_mark_duplicates
----------------------
Runs Picard mark duplicates on the sorted bam file.

metrics
-------
This step runs a series of multiple metrics collection jobs and the output bam
from mark duplicates.

- Picard CollectMultipleMetrics: A collection of picard metrics that runs at the
same time to save on I/O.
    - CollectAlignmentSummaryMetrics,
    - CollectInsertSizeMetrics,
    - QualityScoreDistribution,
    - MeanQualityByCycle,
    - CollectBaseDistributionByCycle
- BVATools DepthOfCoverage: Using the specified `BED Files` in the sample sheet,
calculate the coverage of each target region.
- Picard CalculateHsMetrics: Calculates a set of Hybrid Selection specific
metrics from the BAM file. The bait and interval list is automatically created
from the specicied `BED Files`.

md5
---
Create md5 checksum files for the fastq, bam and bai using the system 'md5sum'
util.

One checksum file is created for each file.

report
------
Generate a JSON file reporting the whole pipeline.
The jobs of this step actually update the JSON report as the pipeline is running

copy
----
Copy the whole processing folder to where they can be serve or loaded into a LIMS

final_notification
------------------
Writes a simple '.done' file when the whole pipeline is over

basecall
--------
Use write_fastq software from MGI to perform the base calling.
Takes the raw .cal files from the sequencer and produces fastq files.
Demultiplexing with MGI splitBarcode while doing the basecalling can
be perform if requested with --splitbarcode-demux


