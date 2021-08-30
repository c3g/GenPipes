[TOC]


Illumina Run Processing Pipeline
================================

The standard MUGQIC Illumina Run Processing pipeline uses the Illumina bcl2fastq
software to convert and demultiplex the base call files to fastq files. The
pipeline runs some QCs on the raw data, on the fastq and on the alignment.

Sample Sheets
-------------

The pipeline uses two input sample sheets. The first one is the standard Casava
sheet, a csv file having the following columns (please refer to the Illumina
Casava user guide):

- `SampleID`
- `FCID`
- `SampleRef`
- `Index`
- `Description`
- `Control`
- `Recipe`
- `Operator`
- `SampleProject`

Example:

    FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject
    H84WNADXX,1,sample1_MPS0001,,TAAGGCGA-AGAGTAGA,,N,,,nanuq
    H84WNADXX,1,sample47_MPS0047,,GTAGAGGA-CTAAGCCT,,N,,,nanuq


The second sample sheet is called the Nanuq run sheet. It's a csv file with the
following minimal set of mandatory columns (the column order in the file doesn't
matter)

- `ProcessingSheetId` Must be the same as the `SampleID` from the Casava Sheet.
- `Name` The sample name put in RG headers of bam files and on filename on disk.
- `Run` The run number.
- `Region` The lane number.
- `Library Barcode` The library barcode put in .bam's RG headers and on disk
- `Library Source` The type of library. If this value contains `RNA` or `cDNA`,
`STAR` will be used to make the aligmnent, otherwise, `bwa_mem` will be used
- `Library Type` Used to determine is the sample is from cDNA/RNA when the
`Library Source` is `Library`
- `BED Files` The name of the BED file containing the genomic targets. This is
the `filename` parameter passed to the `fetch_bed_file_command`
- `Genomic Database` The reference used to make the alignment and calculate aligments metrics

Example:

    Name,Genomic Database,Library Barcode,Library Source,Library Type,Run,Region,BED Files,ProcessingSheetId
    sample1,Rattus_norvegicus:Rnor_5.0,MPS0001,RNA,Nextera XT,1419,1,toto.bed,sample1_MPS0001
    sample47,,MPS1047,Library,Nextera XT,1419,2,toto.bed,sample47_MPS1047


Usage
-----
```
#!text

usage: illumina_run_processing.py [-h] [--help] [-c CONFIG [CONFIG ...]]
                                  [-s STEPS] [-o OUTPUT_DIR]
                                  [-j {pbs,batch,daemon,slurm}] [-f]
                                  [--no-json] [--report] [--clean]
                                  [-l {debug,info,warning,error,critical}]
                                  [--sanity-check]
                                  [--container {wrapper, singularity} <IMAGE PATH>]
                                  [--genpipes_file GENPIPES_FILE] [-d RUN_DIR]
                                  [--lane LANE_NUMBER] [-r READSETS]
                                  [-i CASAVA_SHEET_FILE] [-x FIRST_INDEX]
                                  [-y LAST_INDEX] [-m NUMBER_OF_MISMATCHES]
                                  [-w] [-v]

Version: 3.6.0

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
                        Run inside a container providing a validsingularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -d RUN_DIR, --run RUN_DIR
                        run directory
  --lane LANE_NUMBER    lane number
  -r READSETS, --readsets READSETS
                        nanuq readset file. The default file is
                        'run.nanuq.csv' in the output folder. Will be
                        automatically downloaded if not present.
  -i CASAVA_SHEET_FILE  illumina casava sheet. The default file is
                        'SampleSheet.nanuq.csv' in the output folder. Will be
                        automatically downloaded if not present
  -x FIRST_INDEX        first index base to use for demultiplexing
                        (inclusive). The index from the sample sheet will be
                        adjusted according to that value.
  -y LAST_INDEX         last index base to use for demultiplexing (inclusive)
  -m NUMBER_OF_MISMATCHES
                        number of index mistmaches allowed for demultiplexing
                        (default 1). Barcode collisions are always checked.
  -w, --force-download  force the download of the samples sheets (default:
                        false)
  -v, --version         show the version information and exit

Steps:
------
1- index
2- fastq
3- align
4- picard_mark_duplicates
5- metrics
6- blast
7- qc_graphs
8- md5
9- copy
10- end_copy_notification

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
Launch fastq generation from Illumina raw data using BCL2FASTQ conversion
software.

The index base mask is calculated according to the sample and run configuration;
and also according the mask parameters received (first/last index bases). The
Casava sample sheet is generated with this mask. The default number of
mismatches allowed in the index sequence is 1 and can be overrided with an
command line argument. No demultiplexing occurs when there is only one sample in
the lane.

An optional notification command can be launched to notify the start of the
fastq generation with the calculated mask.

align
-----
Align the reads from the fastq file, sort the resulting .bam and create an index
of that .bam.

An basic aligment is performed on a sample when the `SampleRef` field of the
Illumina sample sheet match one of the regexp in the configuration file and the
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

blast
-----
Run blast on a subsample of the reads of each sample to find the 20 most
frequent hits.

The `runBlast.sh` tool from MUGQIC Tools is used. The number of reads to
subsample can be configured by sample or for the whole lane. The output will be
in the `Blast_sample` folder, under the Unaligned folder.

qc_graphs
---------
Generate some QC Graphics and a summary XML file for each sample using 
[BVATools](https://bitbucket.org/mugqic/bvatools/).

Files are created in a 'qc' subfolder of the fastq directory. Examples of
output graphic:

- Per cycle qualities, sequence content and sequence length;
- Known sequences (adaptors);
- Abundant Duplicates;

md5
---
Create md5 checksum files for the fastq, bam and bai using the system 'md5sum'
util.

One checksum file is created for each file.

copy
----
Copy processed files to another place where they can be served or loaded into a
LIMS.

The destination folder and the command used can be set in the configuration
file.

An optional notification can be sent before the copy. The command used is in the configuration file.

end_copy_notification
---------------------
Send an optional notification to notify that the copy is finished.

The command used is in the configuration file. This step is skipped when no
command is provided.


