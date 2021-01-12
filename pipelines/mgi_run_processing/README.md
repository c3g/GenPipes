[TOC]


MGI Run Processing Pipeline
================================

The standard MUGQIC MGI Run Processing pipeline uses fastq files produced
by the sequencer, then does demultiplexing. Finally, the
pipeline runs some QCs on the raw data, on the fastq and on the alignment.

Sample Sheets
-------------

The pipeline uses one input sample sheet.
CURRENTLY BASED ON MGI RUN PROCESSING GOOGLE SHEET (https://docs.google.com/spreadsheets/d/1Jk11bQUJdqVg37gfn7ndfk-g9ke96tsjCfAsf3r1xdA)
A csv file having the following columns :

- Sample
- Readset
- Library
- Project
- Project ID
- Protocol
- Index
- Pool ID
- Run ID
- Flowcell ID
- Lane
- Sequencer
- Sequencer ID

Example:
    Sample,Readset,Library,Project,Project ID,Protocol,Index,PoolID,RunID,FlowcellID,Lane,Sequencer,SequencerID
    LSPQ_Viral_Culture_dil_10-1_10cycles,LSPQ_Viral_Culture_dil_10-1_10cycles_PROD_000034-A01,PROD_000034-A01,LSPQ,,CleanPlex_MGI,1,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
    LSPQ_Viral_Culture_dil_10-2_10cycles,LSPQ_Viral_Culture_dil_10-2_10cycles_PROD_000034-B01,PROD_000034-B01,LSPQ,,CleanPlex_MGI,2,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
    LSPQ_Viral_Culture_dil_10-3_10cycles,LSPQ_Viral_Culture_dil_10-3_10cycles_PROD_000034-C01,PROD_000034-C01,LSPQ,,CleanPlex_MGI,3,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
    LSPQ_Nasal_Swab_Neg_ctl_10cycles,LSPQ_Nasal_Swab_Neg_ctl_10cycles_PROD_000034-A02,PROD_000034-A02,LSPQ,,CleanPlex_MGI,25,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
    LSPQ_Viral_Culture_dil_10-4_13cycles,LSPQ_Viral_Culture_dil_10-4_13cycles_PROD_000034-A03,PROD_000034-A03,LSPQ,,CleanPlex_MGI,28,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
    LSPQ_Viral_Culture_dil_10-5_13cycles,LSPQ_Viral_Culture_dil_10-5_13cycles_PROD_000034-B03,PROD_000034-B03,LSPQ,,CleanPlex_MGI,29,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
    L00241026_dil_10-1_13cycles,L00241026_dil_10-1_13cycles_PROD_000034-E03,PROD_000034-E03,LSPQ,,CleanPlex_MGI,33,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
    L00241026_dil_10-2_13cycles,L00241026_dil_10-2_13cycles_PROD_000034-F03,PROD_000034-F03,LSPQ,,CleanPlex_MGI,34,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
    Arctic_RT_reaction_13cycles,Arctic_RT_reaction_13cycles_PROD_000034-B04,PROD_000034-B04,LSPQ,,CleanPlex_MGI,4,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01



Usage
-----
```
#!text

usage: mgi_run_processing.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                             [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}]
                             [-f] [--no-json] [--report] [--clean]
                             [-l {debug,info,warning,error,critical}]
                             [--sanity-check]
                             [--container {wrapper, singularity} <IMAGE PATH>]
                             [-r READSETS] [-d RUN_DIR] [--run-id RUN_ID]
                             [--flowcell-id FLOWCELL_ID]
                             [--raw-fastq-prefix RAW_FASTQ_PREFIX]
                             [--lane LANE_NUMBER] [--demux-fastq]
                             [-x FIRST_INDEX] [-y LAST_INDEX]
                             [-m NUMBER_OF_MISMATCHES]
                             [--allow-barcode-collision] [-v]

Version: 3.1.6-beta

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
  -r READSETS, --readsets READSETS
                        Sample sheet for the MGI run to process (mandatory)
  -d RUN_DIR, --run RUN_DIR
                        Run directory (mandatory)
  --run-id RUN_ID       Run ID. Default is parsed from the run folder
  --flowcell-id FLOWCELL_ID
                        Flowcell ID. Default is parsed from the run folder
  --raw-fastq-prefix RAW_FASTQ_PREFIX
                        Prefix used to search for the raw fastq from the
                        sequencer. Default <FLOWCELL_ID>_<RUN_ID>
  --lane LANE_NUMBER    Lane number (to only process the given lane)
  --demux-fastq         Fastq files given by the sequencer are already
                        demultiplexed : NO DEMULTIPLEXING will be performed by
                        the pipeline
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
1- index
2- fastq
3- qc_graphs
4- fastqc
5- blast
6- align
7- picard_mark_duplicates
8- metrics
9- report
10- copy
11- final_notification

```
index
-----
First copy all the files of the lane from the sequencer deposit folder to the processing folder,
into "raw_fastq".
Then, if demultiplexing already done on sequencer, formerly the case when MGI adapters were used, then
    rename the fastq files and move them from "raw_fastq" to "Unaligned.LANE" folder,
    --> pipeline will SKIP demultiplexing (i.e. next step, fastq).
Else, (i.e. demultiplexing still remains to be done) nothing remains to do here : 
    everything is in place for the demultiplexing to happen in the next step
*TO DO* - in both cases, retrieve index stats from the sequencer output files to build a proper index report

fastq
-----
*** In the future, may generate the fastq files from the raw CAL files. ***
Perform demultplexing of the reads with fgbio DemuxFastqs

qc_graphs
---------
Generate some QC Graphics and a summary XML file for each sample using 
[BVATools](https://bitbucket.org/mugqic/bvatools/).

Files are created in a 'qc' subfolder of the fastq directory. Examples of
output graphic:

- Per cycle qualities, sequence content and sequence length;
- Known sequences (adaptors);
- Abundant Duplicates;

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

report
------
Generate a JSON file reporting the whole pipeline

copy
----
Copy the whole processing foler to where they can be serve or loaded into a LIMS

final_notification
------------------
Writes a simple '.done' file when all pipeline is done processing


